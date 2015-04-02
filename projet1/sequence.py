"""
This module contain an abstract data type: Sequence
which represent a sequence of amino acid used in bioinformatic.
"""


class Sequence(list):
    """
    This class is inherited from list and represent a list
    of amino acid which can be used to calculate the alignment
    of two sequences.
    """

    def __repr__(self):
        return "".join(self)

    def __str__(self):
        return "".join(self)

    @staticmethod
    def fromFile(files):
        """
        This method take a path to a file containing a amino acid sequence
        and generate Sequences from this file (generator).
        """
        data = ""
        for line in files:
            if line[0] == '>':
                break
        for line in files:
            if len(line) > 0:
                if line[0] != '>':
                    data += line.strip("\n").replace(" ", "")
                else:
                    yield Sequence(data)
                    data = ""
