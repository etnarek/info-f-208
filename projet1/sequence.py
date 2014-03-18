class Sequence(list):

    def __repr__(self):
        return "".join(self)

    def __str__(self):
        return "".join(self)

    @staticmethod
    def fromFile(filename):
        with open(filename) as f:
            data = ""
            line = f.readline()
            if line[0] != '>':
                data += line.strip("\n")
            for line in f:
                if line[0] != '>':
                    data += line.strip("\n")
                else:
                    yield Sequence(data)
                    data = ""
