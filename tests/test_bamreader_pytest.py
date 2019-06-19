from msi.bamreader import BamReader
from msi.bedreader import BedReader
from msi.collect_msi import collect_fasta_sites
from msi.collect_msi import get_sites_from_bam
from msi.collect_msi import compare_possible_msi
from msi.collect_msi import extend_kmers
from msi.kmers import Kmers

import pytest


def test_bam_test_read():
    path = '../data/bed_bam_test/test.bam'
    bam = BamReader(path).file
    bam_iter = bam.fetch('1', 28234090, 28234093)
    list1 = [x for x in bam_iter]
    assert len(list1) == 4


def test_bed_test_read_correct():
    path = '../data/bed_bam_test/test1.bed'
    exons = BedReader(path).exonList
    assert len(exons) == 4

    path = '../data/bed_bam_test/test2.bed'
    exons = BedReader(path).exonList
    assert len(exons) == 4


def test_bed_test_read_incorrect():
    path = '../data/bed_bam_test/test3.bed'
    with pytest.raises(IndexError):
        BedReader(path)


def test_collect_fasta_sites():
    bed_path = '../data/bed_bam_test/hg19.bed'
    exons = BedReader(bed_path).exonList
    fasta_path = '/home/test/Share1/ngs/fa/hg19.fa'

    kmers = Kmers().kmers
    for exon in exons:
        collect_fasta_sites(kmers, fasta_path, exon)


# TODO: split to unit tests
def test_collect_exome_chr7_fasta_sites():
    bed_path = '../data/bed_bam_test/panel_az600_chr7.bed'
    exons = BedReader(bed_path).exonList
    fasta_path = '/home/test/Share1/ngs/fa/hg38.fa'

    tumor = '../data/bed_bam_test/Dev_731_GTL_16_5_Pool1-ready.bam'
    normal = '../data/bed_bam_test/Dev_731_NA12878a_Pool1-ready.bam'

    msis = []
    kmers = Kmers().kmers
    for exon in exons:
        position_to_kmer, reference = collect_fasta_sites(kmers, fasta_path, exon)
        extended_kmers = extend_kmers(position_to_kmer, reference)
        possible_msi_tumor = get_sites_from_bam(tumor, reference, extended_kmers)
        possible_msi_normal = get_sites_from_bam(normal, reference, extended_kmers)
        msi = compare_possible_msi(possible_msi_normal, possible_msi_tumor)
        if len(msi) > 0:
            msis.append(msi)
    print(msis)
