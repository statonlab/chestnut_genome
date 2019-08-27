"""
Created on Aug 27th, 2019
Author: Jiali

This program takes blast results and add the chromosome infomation to each genes based on the gtf/gff files.
"""
import argparse
parser = argparse.ArgumentParser(description="Add chromosome info to genes", usage="%(prog)s -i input_BLAST_file -o output_filename -gff chestnutGFF_file -gtf oakGTF_file")
parser.add_argument("-i", type=str, help="input blast tabular file")
parser.add_argument("-gff", type=str, help="gff annotation file")
parser.add_argument("-gtf", type=str, help="gtf annotation file")
parser.add_argument("-o", type=str, help="output file")
args = parser.parse_args()
BLASTfile = args.i
output = args.o
chestnutGFF = args.gff
oakGTF = args.gtf

def getGenePosition(id, gff):
    with open(gff) as annotation:
        for line in annotation:
            if id in line:
                positions = line.strip("\n").split("\t")
                contig = positions[0]
                start = positions[3]
                end = positions[4]
                return contig+'\t'+start+'\t'+end
def main():
    with open(BLASTfile) as file_in, open(output,"w") as file_out:
        for line in file_in:
            pairs = line.strip("\n").split("\t")
            CmID = pairs[0]
            QlID = pairs[1]
            identity = pairs[3]
            CmPosition = getGenePosition(CmID, chestnutGFF)
            QlPosition = getGenePosition(QlID, oakGTF)
            file_out.write(str(CmID)+'\t'+str(CmPosition)+'\t'+str(QlID)+'\t'+str(QlPosition)+'\t'+str(identity)+'\n')
main()