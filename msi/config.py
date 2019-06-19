# Configuration parameters for thresholds of quality etc.


class Config:
    freq = 0.05
    mapq = 10.0
    baseq = 25
    mincov = 2
    th = 1
    log_level = 'WARNING'

    @staticmethod
    def set_config(args):
        if args.freq:
            Config.freq = args.freq
        if args.mapq:
            Config.mapq = args.mapq
        if args.baseq:
            Config.baseq = args.baseq
        if args.mincov:
            Config.mincov = args.mincov
        if args.th:
            Config.th = args.th
        if args.log:
            Config.log_level = args.log

    def __init__(self, args):
        self.normal = args.normal
        self.fasta = args.fasta
        self.tumor = args.tumor
        self.bed = args.bed
        self.set_config(args)
