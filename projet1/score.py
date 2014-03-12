sortList = [
    'A',
    'R',
    'N',
    'D',
    'C',
    'Q',
    'E',
    'G',
    'H',
    'I',
    'L',
    'K',
    'M',
    'F',
    'P',
    'S',
    'T',
    'W',
    'Y',
    'V',
    'B',
    'Z',
    'X',
    '*']


class Score(dict):

    def __int__(self, *arg, **kw):
        dict.__init__(*arg, **kw)

    def __repr__(self):
        string = "    "
        string += "   ".join(sortList) + "\n"
        for i in sortList:
            string += i
            for j in sortList:
                string += "{0:>4}".format(self[i,j])
            string += "\n"
        return string

    def __str__(self):
        return repr(self)

    @staticmethod
    def fromFile(filename):
        matrice = Score()
        with open(filename) as f:
            for line in f:
                if line[0] == ' ':
                    up = line.split()
                elif line[0] != '#':
                    head = line[0]
                    data = line[1:].split()
                    for i in range(len(data)):
                        matrice[up[i], head] = data[i]
        return matrice
