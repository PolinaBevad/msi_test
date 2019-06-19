from msi.fastareader import FastaReader
from msi.fastareader import Reference
from msi.bamreader import BamReader


def collect_fasta_sites(kmers, fasta_path, exon):
    reference = FastaReader(fasta_path, exon).reference
    ref_sequence = reference.reference_sequence
    position_to_kmers = {}
    k = 0
    while k < len(ref_sequence):
        for kmer in sorted(kmers, reverse=True):
            i = 0
            j = 0
            while k + i < len(ref_sequence) and j < len(kmer) and ref_sequence[k + i].lower() == kmer[j].lower():
                i += 1
                j += 1

            if j == len(kmer):
                position_to_kmers[k] = kmer
                k += len(kmer) - 1
        k += 1

    return position_to_kmers, reference


def extend_kmers(kmers, reference):
    extended_kmers = {}
    for kmer in kmers:
        extension = 5
        kmer_len = len(kmers[kmer])
        extended = reference.reference_sequence[kmer - extension:kmer + kmer_len + extension].lower()
        possible_kmers = []
        for i in range(kmer_len + 1):
            new_kmer = extended[:extension] + extended[5:kmer_len + extension - i] + extended[kmer_len + extension:]
            possible_kmers.append(new_kmer)
        extended_kmers[kmer - extension] = possible_kmers
    return extended_kmers


def get_sites_from_bam(bam_path, reference, extended_kmers):
    bam = BamReader(bam_path).file
    possible_msi = {}
    for position in extended_kmers:
        reads = bam.fetch(reference.chr, reference.start + position, reference.start + position + 50)
        kmers = {}
        for read in reads:
            for kmer in extended_kmers[position]:
                if kmer.lower() in read.query_alignment_sequence.lower():
                    if kmer in kmers:
                        kmers[kmer] += 1
                    else:
                        kmers[kmer] = 1
                    break
        possible_msi[position + reference.start] = kmers

    return possible_msi


def compare_possible_msi(normal, tumor):
    msi = {}
    for pos in tumor:
        for kmer in tumor[pos]:
            if kmer not in normal[pos]:
                msi[pos] = (kmer, tumor[pos][kmer])
    return msi
