
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
        self.len = len(seq)


class Gene:
    def __init__(self, name:str, seq:str):
        '''Represents one fasta record'''
        ## Data ##
        self.name = name # should be extracted from record header  
        self.seq = seq
        self.exon_start, self.exon_end = self._extract_exon_pos() 
        self.seq = self.seq.lower() # once we have exon position, set the whole sequence to lower case
        self.motif_hits:list = [] # list of MotifInstances, instantiated to an empty list since no hits have been found when object is created

    ## Methods ##
    def _extract_exon_pos(self):
        # 0-based, inclusive
        start = -1
        for i, char in enumerate(self.seq):
            if (char.isupper()):
                start = i
                break
        if (start == -1):
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
    
    def add_motif_hit(self, hit:"MotifInstance"):
        self.motif_hits.append(hit)


class MotifInstance:
    def __init__(self, motifType:Motif, gene:Gene, pos:int):
        '''A motif instance/hit detected in a sequence from the fasta file'''
        ## Data ##
        self.motifType = motifType
        self.gene = gene
        self.pos = pos # 0-based, inclusive 


    
        

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
    """
    Parses a fasta to extract each record as a Gene object, containing name (header line) and sequence. 
    Lowercase bases are considered introns and uppercase bases are considered exons.

    :param fasta: explicit path to fasta file (each record should be an intron-exon-intron region)
    :type fasta: str
    """
    genes:set = set() # set of genes found in fasta
    name:str = ""
    seq:str = ""
    with open(fasta, "r") as file:
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                if name != "":
                    gene = Gene(name, seq)
                    genes.add(gene)
                name = line[1:-1] # extract header line as name 
            else:
                seq += line
    return genes


def parse_motifs(motifs:str):
    """
    Parses a text file containing motif sequences to extract each motif as a Gene object, containing motif sequence and regex for this motif.

    :param motifs: explicit path to motifs text file (one motif per line)
    :type motifs: str
    """
    motif_set:set = set()
    seen_seqs:list = []
    with open(motifs, "r") as file:
        for line in file:
            line = line.strip()
            seq = line.lower() # extract motif sequence, make all bases lowercase for consistency
            if seq in seen_seqs:
                continue # skip any duplicate motifs observed
            seen_seqs.append(seq)
            regex = motif_seq_to_regex(seq) # get regex for motif sequence
            motif_set.add(Motif(seq, regex))
    return motif_set

def motif_seq_to_regex(motif:str):
    """
    Takes a motif sequence and returns a regex string representing that motif. Capitalization of sequence is ignored.
    
    :param motif: motif sequence (may contain ambiguous bases)
    :type motif: str
    """
    regex:str = ""
    nucleotides:dict = {
        "a":"a",
        "t":"t",
        "c":"c",
        "g":"g",
        "u":"u",
        "w":{"a","t","u"},
        "s":{"c", "g"},
        "m":{"a", "c"},
        "k":{"g", "t", "u"},
        "r":{"a", "g"},
        "y":{"c", "t", "u"},
        "b":{"c", "g", "t", "u"},
        "d":{"a", "g", "t", "u"},
        "h":{"a", "c", "t", "u"},
        "v":{"a", "c", "g"},
        "n":{"a","g","c", "t", "u"},
    }

    for char in motif:
        if char == "-":
            print(f"At least one gap ('-') found in the motif '{motif}'. Motifs must not contain any gaps.")
            sys.exit(1)
        if char not in nucleotides:
            print(f"Invalid character '{char}' encountered in the motif '{motif}'. Only nucleotides (a,t,c,g,u) or amibiguous base symbols are accepted.")
            sys.exit(1)
        valid_nts = nucleotides[char] # get valid nucleotides for this base
        if len(valid_nts) == 1: # if nt is an unamiguous base 
            regex += char
        else: # if nt is an ambiguous base 
            regex += "[" + "".join(valid_nts) + "]"

    return regex

def find_motif_hits(gene:Gene, motif_set:set[Motif]):
    """
    Finds all motif instances in a gene and adds them to the gene object's motif_hits list.
    
    :param gene: Gene (aka record) to find motif hits in 
    :type gene: Gene
    :param motif_set: Motifs of interest
    :type motif_set: set of Motifs
    """
    
    for motif in motif_set:
        regex:str = motif.regex
        seq:str = motif.seq
        for match in re.finditer(f"(?=({regex}))", seq): # use lookahead to make sure we find overlapping motifs 
            hit:MotifInstance = MotifInstance(motif, gene, match.start())
            gene.add_motif_hit(hit)

def make_plot(genes:list[Gene]):
    """
    Create one plot for all records in the fasta using pycairo.
    
    :param genes: Genes (records) in the fasta, each with a populated motif_hits list 
    :type genes: list[Gene]
    """

    for gene in genes:
        

def Main():
   args = get_args()
   fasta = args.fasta
   motifs = args.motifs
   genes = parse_fasta(fasta) # set of Gene objects extracted from fasta
   motif_set = parse_motifs(motifs) # set of Motif objects extracted from motifs file
   for gene in genes: # populate each gene's motif_hits list with identified motifs 
       find_motif_hits(gene, motif_set)




if __name__ == "__main__":
   Main()
    
