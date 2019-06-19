class Kmers:

    def __init__(self):
        self.kmers = []
        mono_min_length = 5
        mono_max_length = 10

        di_min_length = 3
        di_max_length = 6

        bases = ['A', 'C', 'G', 'T']

        for base in bases:
            for i in range(mono_min_length, mono_max_length + 1):
                self.kmers.append("".join(base) * i)

        for base1 in bases:
            for base2 in bases:
                if base1 == base2:
                    continue
                for i in range(di_min_length, di_max_length + 1):
                    self.kmers.append((base1 + base2) * i)

    def __str__(self):
        return str(self.kmers)
