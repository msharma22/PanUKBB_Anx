#!/usr/bin/python3
'''This python script takes one gwas summary statistics file
as input, retrieves dbSNP153 rsids, prepends gwas summary statistics
 as output'''

import gzip
import re
import sys
import argparse
import os

def check_arg(args=None):
    parser = argparse.ArgumentParser(description='Script to retrieve rsids')
    parser.add_argument('-i', '--input',
                        help='path to GWAS summary statistics file',
                        required='True'
                        )
    parser.add_argument('-o', '--output',
                        help='output file prefix',
                        required='True'
                        )
    parser.add_argument('-c', '--chr',
                        help='chromosome column name',
                        type=str,
                        required='True'
                        )
    parser.add_argument('-p', '--pos',
                        help='position column name',
                        type=str,
                        required='True'
                        )
    parser.add_argument('-b', '--build',
                        help='genome build, 37 or 38, default 37',
                        type=int,
                        default=37)

    return parser.parse_args(args)

#retrieve command line arguments

args = check_arg(sys.argv[1:])
infile = args.input
outfile = args.output
chr_col = args.chr
pos_col = args.pos
build = args.build

#function to retrieve index of columns containing chr and pos info
#this will be used to pull the chr and pos of each line
def get_index(header, chr, pos):
    for i,head in enumerate(header):
        if head.decode("utf-8")  == chr:
            chr_i = i
        elif head.decode("utf-8") == pos:
            pos_i = i
    return chr_i, pos_i

outgwas = gzip.open(outfile + ".txt.gz","wb")

with gzip.open(infile) as f:
    #just read the first line to retrieve index of 
    #chr and pos in header row
    firstline = f.readline().rstrip().split()
    chr_index, pos_index = get_index(firstline, chr_col, pos_col)
    #write firstline to fine
    outgwas.write(b"rsid\t" + b"\t".join(firstline) + b"\n") #get bytes with b"
    for line in f:
        arr = line.strip().split()
        c = arr[chr_index].decode("utf-8")
        p = arr[pos_index].decode("utf-8")
        #subtract 1 from position
        p1 = str(int(p) - 1)
        #retrieve dbSNP info
        os.system("bigBedToBed /home/isabelle/dbSnp155.bb.1 -chrom=chr" + c + " -start=" + p1 + " -end=" + p + " tmp.txt")
        dbsnplist = open("tmp.txt").read().split()
        if len(dbsnplist) > 0:
            rsid = dbsnplist[3].encode('utf-8') #encode string as bytes
        else:
            rsid = "NA".encode('utf-8') #encode string as bytes
        outgwas.write(rsid + b"\t" + line) #get tab as bytes with b"\t"

outgwas.close()
