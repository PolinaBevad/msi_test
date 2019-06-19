import pysam


class FastaReader:
    EXTENSION = 50

    def __init__(self, path, exon):
        self.reference = self.read_fasta(path, exon)

    def read_fasta(self, path, exon):
        fasta_file = pysam.FastaFile(path)
        start = exon.start - FastaReader.EXTENSION if exon.start > FastaReader.EXTENSION else exon.start
        end = exon.end + FastaReader.EXTENSION \
            if fasta_file.get_reference_length(exon.chr) < exon.end + FastaReader.EXTENSION else exon.end
        reference = fasta_file.fetch(exon.chr, start, end)

        return Reference(reference, exon.chr, start, end)


class Reference:
    def __init__(self, reference, chr, start, end):
        self.reference_sequence = reference
        self.chr = chr
        self.start = start
        self.end = end
