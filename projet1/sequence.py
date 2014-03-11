class Sequence:

    def __init__(self, val = ""):
        self.sequence = val
        self.pos = 0

    def __add__(self, val):
        self.sequence.extend(val)

    def __len__(self):
        return len(self.sequence)

    def __getitem__(self, i):
        return self.sequence[i]

    def __setitem__(self, i, val):
        self.sequence[i] = val

    def __iter__(self):
        self.pos = 0
        return self

    def __next__(self):
        if self.pos > len(self):
            return None
        else:
            pos+=1
            return self[pos - 1]

    @staticmethod
    def fromFile(filename):
        with open(filename) as f:
            data = ""
            line = f.readline()
            if line[0] != '>':
                data+=line.replace("\n","")
            for line in f:
                if line[0] != '>':
                    data+=line.replace("\n","")
                else:
                    yield data
                    data = ""
