"""
This module contain an abstract data type: Score
which represent a matrix (BLOSUM or PAM) used in bioinformatic.
"""
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
    """
    This class is inherited from dict and represent a matrix 2*2.
    This matrix serve to calculate score of an amino acid.
    The wright type of matrix are PAM and BLOSUM.
    """

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
    def fromFile(files):
        """
        This method take a path to a file containing a matrix
        BLOSUM or PAM and return this matrix.
        """
        matrix = Score()
        for line in files:
            if len(line) > 0:
                if line[0] == ' ':
                    firstRaw = line.split()
                elif line[0] != '#':
                    head = line[0]
                    data = line[1:].split()
                    for i in range(len(data)):
                        matrix[firstRaw[i], head] = int(data[i])
        return matrix
