from msi.collect_msi import collect_fasta_sites
from msi.kmers import Kmers
from msi.bedreader import BedReader

def main():
    bed_path = '../data/bed_bam_test/panel_az600_chr7.bed'
    fasta_path = '/home/test/Share1/ngs/fa/hg19.fa'

    exons = BedReader(bed_path).exonList
    kmers = Kmers().kmers

    for exon in exons:
        # collect_fasta_sites()
        position_to_kmer = collect_fasta_sites(kmers, fasta_path, exon)
        print(exon, position_to_kmer)
        # take normal msis or control
        # take somatic msis

