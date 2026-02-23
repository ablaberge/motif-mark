
import cairo
import argparse
import re
import sys
import os

##--------------------------Global Variables--------------------------##

COLORS:list = [
    (1.0, 0.0, 0.0, 0.5), # red
    (0.0, 0.5, 1.0, 0.5), # blue
    (0.0, 0.8, 0.0, 0.5), # green
    (0.6, 0.0, 0.8, 0.5), # purple
    (1.0, 0.6, 0.0, 0.5)  # orange
]

##--------------------------Classes--------------------------##

class Motif:
    def __init__(self, seq:str, regex:str, color:str):
        '''A motif type, extracted from motif text file parameter'''
        ## Data ##
        self.seq = seq 
        self.regex = regex
        self.len = len(seq)
        self.color = color


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
                break
        if (end==-1):
            print("Error: Failed to find the end of the exon. Make sure each record in your fasta includes an intron-exon-intron sequence, introns are denoted with lowercase characters, and exons are uppercase.")
            sys.exit(1)
        #-----------DEBUGGING - delete  before submission------------
        print(f"Gene: {self.name}, Seq len: {len(self.seq)}, Exon: {start}-{end}")
        print(f"First 50 chars: '{self.seq[:50]}'")
        print(f"Chars 240-250: '{self.seq[240:250]}'")
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

def parse_fasta(fasta:str) -> list[Gene]:
    """
    Parses a fasta to extract each record as a Gene object, containing name (header line) and sequence. 
    Lowercase bases are considered introns and uppercase bases are considered exons.

    :param fasta: explicit path to fasta file (each record should be an intron-exon-intron region)
    :type fasta: str
    """
    genes:list[Gene] = [] # list of genes found in fasta
    name:str = ""
    seq:str = ""
    with open(fasta, "r") as file:
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                if name != "":
                    gene = Gene(name, seq)
                    genes.append(gene)
                name = line[1:] # extract header line as name 
                seq = ""
            else:
                seq += line
        if name != "":
            gene = Gene(name, seq)
            genes.append(gene)

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
            motif_set.add(Motif(seq, regex, COLORS[len(seen_seqs)-1]))
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
        for match in re.finditer(f"(?=({regex}))", gene.seq): # use lookahead to make sure we find overlapping motifs 
            hit:MotifInstance = MotifInstance(motif, gene, match.start())
            gene.add_motif_hit(hit)


def make_plot(genes:list[Gene], basename:str, max_seq_len:int, motifs:set[Motif]):
    """
    Create one plot for all records in the fasta using pycairo.
    
    :param genes: Genes (records) in the fasta, each with a populated motif_hits list 
    :type genes: list[Gene]
    :param basename: Basename of fasta, used for naming the output png (e.g. Figure_1.fa -> Figure_1.png)
    :type basename: String
    :param max_seq_len: Length of longest record in fasta
    :type max_seq_len: Integer    
    :param motifs: Set of motif types for making the legend
    :type motifs: Set of Motif objects
    """
    # Each base is treated as 1 pixel 
    num_genes:int = len(genes)
    HEIGHT_MAIRGIN = 60
    LEFT_MARGIN = 100
    HEIGHT_PER_RECORD = 120
    SPACE_BW_RECORDS = 20
    WIDTH:int = 1400
    HEIGHT:int = 2*HEIGHT_MAIRGIN + (num_genes * HEIGHT_PER_RECORD) + ((num_genes - 1) * SPACE_BW_RECORDS)
    LEGEND_HEIGHT = 30 + len(motifs) * 25
    LEGEND_WIDTH = 180
    LEGEND_X = WIDTH - LEGEND_WIDTH - 20 
    LEGEND_Y = 80
    LEGEND_PADDING = 15 
    ENTRY_HEIGHT = 15 
    COLOR_BOX_SIZE = 18 
    BOX_TEXT_GAP = 10
    INTRON_THICKNESS = 2
    EXON_HEIGHT = 20
    MOTIF_HEIGHT = 15 

    


    # Start building graph
    surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, WIDTH, HEIGHT)
    ctx = cairo.Context(surface)
    #ctx.scale(WIDTH, HEIGHT)  # Normalizing the canvas

    # Set background to white 
    ctx.set_source_rgb(1,1,1)
    ctx.rectangle(0, 0, WIDTH, HEIGHT)
    ctx.fill()

    # Make legend box
    ctx.set_source_rgb(0,0,0) # Set color to black for legend outline
    ctx.rectangle(LEGEND_X, LEGEND_Y, LEGEND_WIDTH, LEGEND_HEIGHT)
    ctx.stroke()

    # Make each motif entry in legend 
    curr_y =  LEGEND_Y + LEGEND_PADDING
    for motif in motifs:
        # Draw color box
        ctx.set_source_rgba(*motif.color) 
        ctx.rectangle(LEGEND_X+LEGEND_PADDING, curr_y, COLOR_BOX_SIZE, COLOR_BOX_SIZE)
        ctx.fill()
        # Write motif 
        ctx.set_source_rgb(0,0,0) # black text
        ctx.set_font_size(12)
        ctx.select_font_face("Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
        ctx.move_to(LEGEND_X + LEGEND_PADDING + COLOR_BOX_SIZE + BOX_TEXT_GAP, curr_y + (COLOR_BOX_SIZE/2) + 4)
        ctx.show_text(motif.seq.upper())
        # Adjust height for next entry 
        curr_y += ENTRY_HEIGHT
    
    # Make graph for each record 

    gene_index = 0
    for gene in genes:
        # Write title of gene
        gene_y = HEIGHT_MAIRGIN + gene_index * (HEIGHT_PER_RECORD + SPACE_BW_RECORDS)
        ctx.set_source_rgb(0,0,0) # Black text
        ctx.set_font_size(14)
        ctx.select_font_face("Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_BOLD)
        ctx.move_to(LEFT_MARGIN, gene_y - 10) # Write title 10 px above gene
        ctx.show_text(gene.name)
        # Make line to denote intron
        gene_line_y = gene_y + 20 # Place gene graph below gene title
        ctx.set_source_rgb(0,0,0)
        ctx.set_line_width(INTRON_THICKNESS)
        ctx.move_to(LEFT_MARGIN, gene_line_y)
        ctx.line_to(LEFT_MARGIN + len(gene.seq), gene_line_y)
        ctx.stroke()
        # Make box to denote exon
        ctx.set_source_rgb(0,0,0)
        exon_top = gene_line_y - (EXON_HEIGHT/2) # Center exon on line 
        ctx.rectangle(gene.exon_start + LEFT_MARGIN + 1, exon_top, gene.exon_end - gene.exon_start, EXON_HEIGHT)
        ctx.fill()
        # Draw motifs
        for motif_hit in gene.motif_hits:
            # Get length of the motif 
            motif_start = motif_hit.pos 
            motif_width = motif_hit.motifType.len 
            # Set color based on motifType's color 
            ctx.set_source_rgba(*motif_hit.motifType.color)
            # Draw rectangle to represent motif
            motif_y = gene_line_y - MOTIF_HEIGHT
            ctx.rectangle(LEFT_MARGIN + motif_start, motif_y, motif_width, MOTIF_HEIGHT)
            ctx.fill()
    
        gene_index += 1

    surface.write_to_png(f"{basename}.png")

def get_longest_seq(genes:list[Gene]):
    """
    Returns the length of the longest sequence in genes

    :param genes: Genes (records) in the fasta, each with a populated motif_hits list 
    :type genes: list[Gene]
    :return max_seq_len: Length of longest record in fasta 
    :type max_seq_len: Integer
    """
    max_seq_len:int = 0 
    for gene in genes:
        if len(gene.seq) > max_seq_len:
            max_seq_len = len(gene.seq)
    return max_seq_len
        

def Main():
   args = get_args()
   fasta = args.fasta
   motifs = args.motifs
   genes:list[Gene] = parse_fasta(fasta) # set of Gene objects extracted from fasta
   motif_set:set[Motif] = parse_motifs(motifs) # set of Motif objects extracted from motifs file
   for gene in genes: # populate each gene's motif_hits list with identified motifs 
       find_motif_hits(gene, motif_set)
   make_plot(genes, os.path.splitext(fasta)[0], get_longest_seq(genes), motif_set)




if __name__ == "__main__":
   Main()
    
