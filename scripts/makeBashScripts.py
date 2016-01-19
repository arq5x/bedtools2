#!/usr/bin/env python
# encoding: utf-8
"""
makeBashScripts.py
"""

import sys
import os


def main():
    tool_map =  {'annotate': 'annotateBed', 
                 'bamtobed': 'bamToBed',
                 'bamtofastq': 'bamToFastq', 
                 'bed12tobed6': 'bed12ToBed6', 
                 'bedpetobam': 'bedpeToBam', 
                 'bedtobam': 'bedToBam', 
                 'closest': 'closestBed',
                 'cluster': 'clusterBed',
                 'complement': 'complementBed',
                 'coverage': 'coverageBed', 
                 'expand': 'expandCols',
                 'flank': 'flankBed', 
                 'genomecov': 'genomeCoverageBed', 
                 'getfasta': 'fastaFromBed',
                 'groupby': 'groupBy',
                 'igv': 'bedToIgv', 
                 'intersect': 'intersectBed', 
                 'links': 'linksBed', 
                 'map': 'mapBed', 
                 'maskfasta': 'maskFastaFromBed', 
                 'merge': 'mergeBed', 
                 'multicov': 'multiBamCov', 
                 'multiinter': 'multiIntersectBed', 
                 'nuc': 'nucBed',
                 'overlap': 'getOverlap', 
                 'pairtobed': 'pairToBed', 
                 'pairtopair': 'pairToPair', 
                 'random': 'randomBed', 
                 'shift': 'shiftBed',
                 'shuffle': 'shuffleBed', 
                 'slop': 'slopBed', 
                 'sort': 'sortBed', 
                 'subtract': 'subtractBed', 
                 'tag': 'tagBam', 
                 'unionbedg': 'unionBedGraphs', 
                 'window': 'windowBed',
                 'makewindows': 'windowMaker'}

    # create a BASH script for each old tool, mapping to the new CLI command.
    for tool in tool_map:
        new = tool
        old = tool_map[tool]
        
        script = open('bin/'  + old, 'w')
        script.write("#!/bin/sh\n")
        script.write("${0%/*}/bedtools " + new + " \"$@\"\n")
        script.close()

if __name__ == "__main__":
    main()