#!/usr/bin/env python

from msi.collect_msi import collect_fasta_sites
from msi.collect_msi import get_sites_from_bam
from msi.collect_msi import compare_possible_msi
from msi.collect_msi import extend_kmers
from msi.kmers import Kmers
from msi.bedreader import BedReader
from msi.config import Config
import argparse


#TODO: add logger
def main():
    config = parse_config()
    exons = BedReader(config.bed).exonList

    kmers = Kmers().kmers
    msis = []
    for exon in exons:
        position_to_kmer, reference = collect_fasta_sites(kmers, config.fasta, exon)
        extended_kmers = extend_kmers(position_to_kmer, reference)
        possible_msi_tumor = get_sites_from_bam(config.tumor, reference, extended_kmers)
        possible_msi_normal = get_sites_from_bam(config.normal, reference, extended_kmers)
        msi = compare_possible_msi(possible_msi_normal, possible_msi_tumor)
        if len(msi) > 0:
            msis.append(msi)
    print(msis)


def parse_config():
    parser = argparse.ArgumentParser()
    required = parser.add_argument_group("required")
    optional = parser.add_argument_group("optional")
    optional.add_argument('--freq', type=float,
                        help="Minimum difference of frequencies for base to consider site as somatic. "
                             "Default: 5%% (0.05)")
    optional.add_argument('--mapq', type=float, help="Minimum read mapping quality. Default: 10.0")
    optional.add_argument('--baseq', type=int, help="Minimum base quiality. Default: 25")
    optional.add_argument('--mincov', type=int,
                        help="Minimum coverage of base for nucleotide to be considered as somatic. "
                             "Default: 2")
    optional.add_argument('--th', type=int, help="Number of threads for multiprocessing mode. Default: 1")
    required.add_argument('--normal', required=True, type=str, help="Path to normal SAM/BAM/CRAM file.")
    required.add_argument('--tumor', required=True, type=str, help="Path to tumor SAM/BAM/CRAM file.")
    required.add_argument('--bed', required=True, type=str, help="Path to BED file.")
    required.add_argument('--fasta', required=True, type=str, help="Path to FASTA file.")

    args = parser.parse_args()
    return Config(args)


if __name__ == "__main__":
    main()
