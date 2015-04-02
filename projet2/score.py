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

    def __repr__(self, proteinList=None):
        string = "    "
        string += "   ".join(SORT_LIST) + "\n"
        if proteinList == None:
            proteinList = SORT_LIST
        for i in proteinList:
            string += i
            for j in proteinList:
                string += "{0:>4}".format(self.get((i, j), 0))
            string += "\n"
        return string

    def __str__(self, proteinList=None):
        return self.__repr__(proteinList)

    @staticmethod
    def fromFile(filename):
        """
        This method take a path to a file containing a matrix
        BLOSUM or PAM and return this matrix.
        """
        matrix = Score()
        with open(filename) as f:
            for line in f:
                if len(line) > 0:
                    if line[0] == ' ':
                        firstRaw = line.split()
                    elif line[0] != '#':
                        head = line[0]
                        data = line[1:].split()
                        for i in range(len(data)):
                            matrix[firstRaw[i], head] = int(data[i])
        return matrix
