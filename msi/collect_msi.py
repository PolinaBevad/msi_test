from msi.fastareader import FastaReader
from msi.config import Config
from msi.bamreader import BamReader
import logging


def collect_fasta_sites(kmers, fasta_path, exon):
    logging.info("Retrieving extended reference for exon %s:%s-%s.", exon.chr, exon.start, exon.end)
    reference = FastaReader(fasta_path, exon).reference
    ref_sequence = reference.reference_sequence
    logging.info("Extended reference %s:%s-%s.", reference.chr, reference.start, reference.end)
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
    logging.info("Kmers for reference region %s:%s-%s are collected.", reference.chr, reference.start, reference.end)

    return position_to_kmers, reference


def extend_kmers(kmers, reference):
    extended_kmers = {}
    for kmer in kmers:
        extension = 5
        kmer_len = len(kmers[kmer])
        extended = reference.reference_sequence[kmer - extension:kmer + kmer_len + extension].lower()
        possible_kmers = []
        # TODO: add kmers that can be more then initial
        for i in range(kmer_len + 1):
            new_kmer = extended[:extension] + extended[5:kmer_len + extension - i] + extended[kmer_len + extension:]
            possible_kmers.append(new_kmer)
        extended_kmers[kmer - extension] = possible_kmers
    logging.info("Kmers for reference region %s:%s-%s extended.", reference.chr, reference.start, reference.end)

    return extended_kmers


def get_sites_from_bam(bam_path, reference, extended_kmers):
    bam = BamReader(bam_path).file
    possible_msi = {}
    for position in extended_kmers:
        reads = bam.fetch(reference.chr, reference.start + position, reference.start + position + 50)
        kmers = {}
        for read in reads:
            if bad_read(read):
                continue
            for kmer in extended_kmers[position]:
                index = read.query_alignment_sequence.lower().find(kmer.lower())
                if index > 0:
                    if read.query_qualities[index] < Config.baseq:
                        logging.debug("Base %s in read %s was filtered by base quality.", read.query_qualities[index],
                                      read.query_name)
                        continue
                    if kmer in kmers:
                        kmers[kmer] += 1
                    else:
                        kmers[kmer] = 1
                    break
        possible_msi[position + reference.start] = kmers
    logging.info("Possible msi region in bam: %s prepared for reference region %s:%s-%s.",
                 bam_path, reference.chr, reference.start, reference.end)

    return possible_msi


def compare_possible_msi(normal, tumor, reference):
    msi = {}
    for pos in tumor:
        for kmer in tumor[pos]:
            if kmer not in normal[pos] and tumor[pos][kmer] > Config.mincov:
                msi[pos] = (kmer, tumor[pos][kmer])
    logging.info("Normal and tumor BAM msi compared for reference region %s:%s-%s.", reference.chr, reference.start,
                 reference.end)

    return msi


# Read doesn't fit criteria for quality
def bad_read(read):
    if read.mapping_quality < Config.mapq:
        logging.debug("Read %s was filtered by mapping quality.", read.query_name)
        return True
    if read.is_supplementary:
        logging.debug("Read %s was filtered as supplementary.", read.query_name)
        return True

    return False