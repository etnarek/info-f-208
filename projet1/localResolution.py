#! /bin/python3

import sys
from sequence import Sequence
from score import Score

def createMatrix(seq1, seq2, I, E, score):
    matrix = [[[0, 0, 0] for j in seq2] for i in seq1]
    for i in range(1, len(seq1)):
        matrix[i][0] = [0] * 3
    for i in range(1, len(seq2)):
        matrix[0][i] = [0] * 3

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
            S = max(S, V, W, 0)
            matrix[i][j][0] = S
    return matrix

def findMax(F):
    maximum = 0
    pos = (0, 0)
    for i in range(len(F)-1, 0-1, -1):
        for j in range(len(F[i])-1, 0-1, -1):
            if F[i][j][0] > maximum:
                maximum = F[i][j][0]
                pos = (i, j)
    return pos

def aligne(seq1, seq2, F, S, I, E):
    align1 = ""
    align2 = ""
    alignChar = ""
    i, j = findMax(F)
    identity = 0
    globalScore = 0
    while i > 0 and j > 0:
        if F[i][j] == 0:
            break
        if(F[i][j][0] == F[i - 1][j - 1][0] + S[seq1[i], seq2[j]]):
            align1 = seq1[i] + align1
            align2 = seq2[j] + align2
            globalScore += S[seq1[i], seq2[j]]
            if seq1[i] == seq2[j]:
                alignChar = ":" + alignChar
                identity += 1
            elif S[seq1[i], seq2[j]] >= 0:
                alignChar = "." + alignChar
            else:
                alignChar = " " + alignChar
            i -= 1
            j -= 1

        elif F[i][j][0] == F[i-1][j][0] - I:
            align1 = seq1[i] + align1
            align2 = "-" + align2
            alignChar = " " + alignChar
            globalScore -= I
            i -= 1
        elif F[i][j][0] == F[i-1][j][1] - E:
            align1 = seq1[i] + align1
            align2 = "-" + align2
            alignChar = " " + alignChar
            globalScore -= E
            i -= 1

        elif F[i][j][0] == F[i][j-1][0] - I:
            align1 = "-" + align1
            align2 = seq2[j] + align2
            alignChar = " " + alignChar
            globalScore -= I
            j -= 1
        elif F[i][j][0] == F[i][j-1][2] - E:
            align1 = "-" + align1
            align2 = seq2[j] + align2
            alignChar = " " + alignChar
            globalScore -= E
            j -= 1
    print(align1)
    print(alignChar)
    print(align2)
    print(identity)
    print(globalScore)
    return align1, align2, alignChar, identity, globalScore


def main():
    fileName = sys.argv[1]
    score = Score.fromFile(fileName)
    fileName = sys.argv[2]
    seq = list(Sequence.fromFile(fileName))
    seq1 = seq[1]
    seq2 = seq[2]

    matrix = createMatrix(seq1, seq2, 14, 4, score)
    for i in matrix:
        for j in i:
            print("{0:<5}".format(j[0]), end="")
        print()
    aligne(seq1, seq2, matrix, score, 14, 4)

if __name__ == '__main__':
    main()
