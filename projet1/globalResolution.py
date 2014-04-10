#! /bin/python3

import sys
from sequence import Sequence
from score import Score

def createF(seq1, seq2, I, E, score):
    matrix = [[[0, 0, 0] for j in seq2] for i in seq1]
    for i in range(1, len(seq1)):
        S = matrix[i-1][0][0] - I
        V = matrix[i-1][0][1] - E
        V = max(S, V)
        matrix[i][0] = [V] * 3
    for i in range(1, len(seq2)):
        S = matrix[i-1][0][0] - I
        W = matrix[i-1][0][2] - E
        W = max(S, W)
        matrix[i][0] = [W] * 3

    for i in range(1, len(seq1)):
        for j in range(1, len(seq2)):
            S = matrix[i-1][j][0] - I
            V = matrix[i-1][j][1] - E
            V = max(S, V)
            matrix[i][j][1] = V

            S = matrix[i][j-1][0] - I
            W = matrix[i][j-1][2] - E
            W = max(S, W)
            matrix[i][j][2] = W

            S = matrix[i-1][j-1][0] + score[seq1[i], seq2[j]]
            S = max(S, V, W)
            matrix[i][j][0] = S
    return matrix

def aligne(seq1, seq2, F):
    pass

if __name__ == '__main__':
    fileName = sys.argv[1]
    score = Score.fromFile(fileName)
    fileName = sys.argv[2]
    seq = list(Sequence.fromFile(fileName))
    seq1 = seq[0]
    seq2 = seq[1]

    print(createF(seq1, seq2, 6, 3, score))
