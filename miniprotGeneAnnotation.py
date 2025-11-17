#! /usr/bin/python3
# -*- coding:UTF-8 -*-
"""
===================================================
@author:朱成龙
@desc:fast homology gene annotation based on miniprot
@date:2023-05-23 20:54:55
@Email: zhuchenglong96@gmail.com
===================================================
"""

import sys
import argparse
import re
import os
import pybedtools
import igraph as ig
from multiprocessing import Process


def get_args():
    # if no args input, print help information
    if len(sys.argv) == 1:
        sys.argv.append("-h")
    parser = argparse.ArgumentParser(description="fast homology gene annotation based on miniprot")
    parser.add_argument("-g", "--genome", help="genome file", required=True)
    parser.add_argument("-p", "--pep", help="pep file", required=True)
    parser.add_argument("-t", "--threads", help="threads [20]", default=20)
    parser.add_argument("-o", "--output", help="output file", required=True)
    parser.add_argument("-P", "--prefix", help="prefix [MP]", default="MP")
    args = parser.parse_args()
    return args

def groupGff3(cdsBed_lines):
    groupGeneList = []
    cdsLines = cdsBed_lines.strip().split("\n")
    geneSplit = re.compile(r':::')
    g = ig.Graph()
    nameList = []
    edgeList = []
    for line in cdsLines:
        if len(line) == 1:
            continue
        lines = line.strip().split("\t")
        geneIds = geneSplit.split(lines[4])
        geneIds = list(set(geneIds))
        # delete ''
        if '' in geneIds:
            geneIds.remove('')
        for i in range(len(geneIds)):
            nameList.append(geneIds[i])
            for j in range(i+1, len(geneIds)):

                edgeList.append((geneIds[i], geneIds[j]))

    nameList = list(set(nameList))
    g.add_vertices(nameList)
    g.add_edges(edgeList)

    groupGene = g.connected_components()
    for group in groupGene:
        # get group gene id
        group = [g.vs[i]['name'] for i in group]
        groupGeneList.append(group)

    singleGeneDict = {}
    for group in groupGeneList:
        for id in group:
            singleGeneDict[id] = tuple(group)

    return singleGeneDict


def readMiniprotOut(miniprotOut):
    CDS_dict = {}
    info_dict = {}
    idSearch = re.compile(r'Parent=([^;]+)')
    mRNA_search = re.compile(r'ID=([^;]+);\S+;Positive=([^;]+);\S*Target=(\S+)')
    Frameshift_search = re.compile(r'Frameshift')
    pre_mature_search = re.compile(r'(\*+)')
    cds_lines = ''
    pepSeq = ''
    pepLength = 0 
    pepStop = 0
    pepSeqDict = {}
    with open(miniprotOut, "r") as f:
        for line in f:
            lines = line.strip().split('\t')
            if line.startswith("##PAF"):
                pepLength = int(lines[4]) - int(lines[3]) + 1
                continue
            elif line.startswith("##STA"):
                pepSeq = lines[1]
                pepSeq = pepSeq[:-2]
                test = pre_mature_search.search(pepSeq)
                if test or pepLength < 30:
                    pepStop = 1
                else:
                    pepStop = 0
                continue
            elif line.startswith("#"):
                continue
            if pepStop == 1:
                continue
            if len(lines) == 1:
                pepStop = 1
                continue
            if lines[2] == 'mRNA' and pepStop == 0:
                if Frameshift_search.search(lines[8]):
                    pepStop = 1
                    continue
                id_identity = mRNA_search.search(lines[8])
                geneId = id_identity.group(1)
                localScore = float(id_identity.group(2)) * pepLength 
                parent = id_identity.group(3)
                info_dict.setdefault(geneId, {})
                info_dict[geneId].setdefault("score", localScore)
                info_dict[geneId].setdefault("parent", parent)
            elif lines[2] == 'CDS' and pepStop == 0:
                geneId = idSearch.search(lines[8])
                lines[8] = geneId.group(1)
                geneId = lines[8]
                lines[3] = int(lines[3])
                lines[4] = int(lines[4])
                line2 = '\t'.join([str(i) for i in lines])
                cds_lines += f'{line2}\n'
                lines[8] = f'Parent={geneId};'
                line2 = '\t'.join([str(i) for i in lines])
                CDS_dict.setdefault(geneId, '')
                CDS_dict[geneId] += f'{line2}\n'
                # use bedtools getfasta to get pep sequence
                pepSeqDict.setdefault(geneId, '')


    cdsBed = pybedtools.BedTool(cds_lines, from_string=True)
    cdsBed = cdsBed.sort().merge(c=[7, 9], o=['collapse', 'collapse'], delim=':::')
    cdsBed_lines = str(cdsBed)
    print("read miniprot file success!")
    return cdsBed_lines, CDS_dict, info_dict

def findBestGene(cdsBed_lines, CDS_dict, info_dict):
    # print(cdsBed_lines + '#')
    singleGeneDict = groupGff3(cdsBed_lines)
    bestGene = {}
    for id in singleGeneDict.keys():
        score = info_dict[id]['score']
        parent = info_dict[id]['parent']
        group = singleGeneDict[id]
        bestGene.setdefault(group, [])
        bestGene[group].append((id, score, parent))
    
    newCDSlines = ''

    for group in (bestGene.keys()):
        bestGene[group].sort(key=lambda x:x[1], reverse=True)
        bestGeneId = bestGene[group][0][0]
        bestGene[group] = bestGeneId
        newCDSlines += f'{CDS_dict[bestGeneId]}'
    
    newCDSbed = pybedtools.BedTool(newCDSlines, from_string=True)
    originalCDSbed = pybedtools.BedTool(cdsBed_lines, from_string=True)
    overlapCDS = originalCDSbed.intersect(newCDSbed, wa=True)
    overlapCDS_lines = str(overlapCDS)
    
    subtractCDS = originalCDSbed.subtract(newCDSbed, A=True)
    subtractCDS_lines = str(subtractCDS)
    deleteGene = []

    geneSplit = re.compile(r':::')
    for line in overlapCDS_lines.strip().split("\n"):
        lines = line.strip().split("\t")
        if lines == ['']:
            continue
        geneId = geneSplit.split(lines[4])
        deleteGene += geneId
    
    noOverlapCDS = []
    for line in subtractCDS_lines.strip().split("\n"):
        lines = line.strip().split("\t")
        if lines == ['']:
            subtractCDS_lines = ''
            break
        geneId = geneSplit.split(lines[4])
        geneId = list(set(geneId))
        if '' in geneId:
            geneId.remove('')
        noOverlapCDS += geneId

    needRemove = set(deleteGene) & set(noOverlapCDS)
    if '' in needRemove:
        needRemove.remove('')
    for id in needRemove:
        subtractCDS_lines = re.sub(id, ':::', subtractCDS_lines)

    return newCDSlines, subtractCDS_lines

def getfinalBestGene(cdsBed_lines, CDS_dict, info_dict, output):

    outputCDS = ''
    while True:
        newCDSlines, cdsBed_lines = findBestGene(cdsBed_lines, CDS_dict, info_dict)
        outputCDS += newCDSlines
        if newCDSlines == '' or cdsBed_lines == '':
            break

    with open(f'{output}', 'w') as f:
        f.write(outputCDS)


def miniprot(genome, pep, threads, output, prefix):
    print('miniprot start!')
    if not os.path.exists(f'{genome}.mpi'):
        cmd = f'/data01/zhuchenglong/software/miniprot/miniprot-0.15_x64-linux/miniprot -t {threads} -d {genome}.mpi {genome}'
        P = Process(target=os.system, args=(cmd,))
        P.start()
        P.join()
    else:
        print(f'\033[31m{genome}.mpi exists! skip miniprot index!\033[0m')
        
    if not os.path.exists(f'{output}.miniprotOut'):
        cmd = f'/data01/zhuchenglong/software/miniprot/miniprot-0.15_x64-linux/miniprot --outc=0.75 --outn=2 -N2 -t {threads} --gff --trans -P {prefix} {genome}.mpi {pep} > {output}.miniprotOut'
        P = Process(target=os.system, args=(cmd,))
        P.start()
        P.join()
    else:
        print(f'\033[31m{output}.miniprotOut exists! skip miniprot!\033[0m')
    
    print('miniprot end!')
            
def main():
    args = get_args()
    genome = args.genome
    pep = args.pep
    threads = args.threads
    output = args.output
    prefix = args.prefix
    miniprotOut = f'{output}.miniprotOut'
    miniprot(genome, pep, threads, output, prefix)
    cdsBed_lines, CDS_dict, info_dict = readMiniprotOut(miniprotOut)
    getfinalBestGene(cdsBed_lines, CDS_dict, info_dict, output)

if __name__ == "__main__":
    main()
