
import cairo
import argparse
import math
import re
import sys

##--------------------------Classes--------------------------##

class Motif:
    def __init__(self, seq:str, regex:str):
        '''A motif type, extracted from motif text file parameter'''
        ## Data ##
        self.seq = seq 
        self.regex = regex


    ## Methods ##
    def get_regex(self):
        return self.regex

class Gene:
    def __init__(self, name:str, seq:str, motifHits:list):
        '''Represents one fasta record'''
        ## Data ##
        self.name = name # should be extracted from record header  
        self.seq = seq
        self.exonStart, self.exonEnd = self._extract_exon_pos() 
        self.motifHits = [] # list of MotifInstances

    ## Methods ##
    def _extract_exon_pos(self):
        # 0-based, inclusive
        start = -1
        for i, char in enumerate(self.seq):
            if (char.isupper()):
                start = i
                break
        if (start==-1):
            print("Error: Failed to find the start of the exon. Make sure each record in your fasta includes an intron-exon-intron sequence, introns are denoted with lowercase characters, and exons are uppercase.")
            sys.exit(1)
        # 0-based, exclusive
        end = -1 
        for i in range(start, len(self.seq)):
            if(self.seq[i].islower()):
                end = i
        if (end==-1):
            print("Error: Failed to find the end of the exon. Make sure each record in your fasta includes an intron-exon-intron sequence, introns are denoted with lowercase characters, and exons are uppercase.")
            sys.exit(1)
        return start, end
    
    def add_motif_hit(self, hit:Motif):
        self.motifHits.append(hit)




class MotifInstance:
    def __init__(self, motifType:Motif, gene:Gene, start_pos:int, end_pos:int):
        '''A motif instance/hit detected in a sequence from the fasta file'''
        ## Data ##
        self.motifTyepe = motifType
        self.gene = gene
        self.start_pos = start_pos # 0-based, inclusive 
        self.end_pos = end_pos # 0-based, exclusive

    
        

##--------------------------Functions--------------------------##
def get_args():
    parser = argparse.ArgumentParser(
        description="Program to visualize motifs in DNA/RNA sequence given a fasta with sequences and motifs text file.")
    parser.add_argument("-f", "--fasta",
                        help="Absolute path to fasta file. Each sequence should contain an exon and the concurrent intron on either side of the given exon.", type=str, required=True)
    parser.add_argument("-m", "--motifs",
                        help="Absolute path to a text file containing motifs of interest. Each motif should be listed on its own line.", type=str, required=True)
    return parser.parse_args()

def parse_fasta(fasta:str):
    # To-Do: Add functionality here
    return

def parse_motifs(motifs:str):
    # To-Do: Add functionality here
    return

def Main():
   args = get_args()
   fasta = args.fasta
   motifs = args.motifs



if __name__ == "__main__":
   Main()
    
