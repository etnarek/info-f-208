SORT_LIST = [
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

    def __repr__(self):
        string = "    "
        string += "   ".join(SORT_LIST) + "\n"
        for i in SORT_LIST:
            string += i
            for j in SORT_LIST:
                string += "{0:>4}".format(self[i, j])
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
                    firstRaw = line.split()
                elif line[0] != '#':
                    head = line[0]
                    data = line[1:].split()
                    for i in range(len(data)):
                        matrice[firstRaw[i], head] = data[i]
        return matrice
