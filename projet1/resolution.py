#!/usr/bin/python3

from sequence import Sequence
from score import Score

def crateMatrix(seq1, seq2, openingPenality, extendingPenality, scoreMatrix, local=False):
    """ Create the matrix used to align 2 sequences given.

    parameters: seq1, seq2 two sequences to be aligned
    openingPenality the penality of creating a new gap
    extendingPenality the penality of extending a gap
    scoreMatrix the pam or blosum matrix used to aligne 2 sequences
    local (default faulse) if we want a local alignement

    retrun: the matrix used to make the alignement
    """
    minimal = 0 if local else float("-inf")
    matrix = [[[0, 0, 0, ""] for j in range(len(seq2) + 1)] for i in range(len(seq1) + 1)]
    matrix[0][0] = [0, 0, 0, "E"]
    for i in range(1, len(seq1) + 1):
        score = max(-openingPenality - (i-1)*extendingPenality, minimal)
        matrix[i][0] = [score, score, minimal, "U"]
    for j in range(1, len(seq2) + 1):
        score = max(-openingPenality - (j-1)*extendingPenality, minimal)
        matrix[0][j] = [score, minimal, score, "L"]
    for i in range(1, len(seq1)+1):
        for j in range(1, len(seq2)+1):
            S = matrix[i-1][j][0] - openingPenality
            V = matrix[i-1][j][1] - extendingPenality
            V = max(S, V, minimal)
            matrix[i][j][1] = V

            S = matrix[i][j-1][0] - openingPenality
            W = matrix[i][j-1][2] - extendingPenality
            W = max(S, W, minimal)
            matrix[i][j][2] = W

            S = matrix[i-1][j-1][0] + scoreMatrix[seq1[i - 1], seq2[j - 1]]
            if S >= V and S >= W:
                matrix[i][j][3] = "D"
            elif V >= W:
                matrix[i][j][3] = "L"
            else:
                matrix[i][j][3] = "U"

            S = max(S, V, W, minimal)
            matrix[i][j][0] = S
    return matrix

