#!/usr/bin/python3

import sys
from sequence import Sequence
from score import Score
import argparse

def createMatrix(seq1, seq2, openingPenality, extendingPenality, scoreMatrix, local=False):
    """ Create the matrix used to align 2 sequences given.

    Parameters: seq1, seq2 two sequences to be aligned
    openingPenality the penalty of creating a new gap
    extendingPenality the penalty of extending a gap
    scoreMatrix the pam or blosum matrix used to align 2 sequences
    local (default false) if we want a local alignment

    Return: the matrix used to make the alignment
    """

    minimal = 0 if local else float("-inf")
    matrix = [[[0, 0, 0, ""] for j in range(len(seq2) + 1)] for i in range(len(seq1) + 1)]
    matrix[0][0] = [0, 0, 0, "E"]
    for i in range(1, len(seq1) + 1):
        score = max(-openingPenality - (i)*extendingPenality, minimal)
        matrix[i][0] = [score, score, minimal, "U"]
    for j in range(1, len(seq2) + 1):
        score = max(-openingPenality - (j)*extendingPenality, minimal)
        matrix[0][j] = [score, minimal, score, "L"]
    for i in range(1, len(seq1)+1):
        for j in range(1, len(seq2)+1):
            S = matrix[i-1][j][0] - (openingPenality + extendingPenality)
            V = matrix[i-1][j][1] - extendingPenality
            V = max(S, V, minimal)
            matrix[i][j][1] = V

            S = matrix[i][j-1][0] - (openingPenality + extendingPenality)
            W = matrix[i][j-1][2] - extendingPenality
            W = max(S, W, minimal)
            matrix[i][j][2] = W

            S = matrix[i-1][j-1][0] + scoreMatrix[seq1[i-1], seq2[j-1]]
            if S == max(S, V, W, minimal):
                matrix[i][j][3] += "D"
            if V == max(S, V, W, minimal):
                matrix[i][j][3] += "L"
            if W == max(S, V, W, minimal):
                matrix[i][j][3] += "U"
            if minimal == max(S, V, W, minimal):
                matrix[i][j][3] += "E"

            S = max(S, V, W, minimal)
            matrix[i][j][0] = S
    return matrix

def bactrackAligne(seq1, seq2, matrix, scoreMatrix, i, j):
    """
    Return all possibles paths from the position i, j until the first E following the U, L and D directions.
    """
    paths = []
    if "E" in matrix[i][j][3]:
        return [["", "", ""]]
    if "D" in matrix[i][j][3]:
        s1 = seq1[i-1]
        s2 = seq2[j-1]
        char = ""
        if s1 == s2:
            char = ":"
        elif scoreMatrix[s1, s2] >= 0:
            char = "."
        else:
            char = " "
        pats = bactrackAligne(seq1, seq2, matrix, scoreMatrix, i-1, j-1)
        for pat in pats:
            pat[0] = pat[0] + s1
            pat[1] = pat[1] + s2
            pat[2] = pat[2] + char
            paths.append(pat)
    if "U" in matrix[i][j][3]:
        s1 = "-"
        s2 = seq2[j-1]
        char = " "
        pats = bactrackAligne(seq1, seq2, matrix, scoreMatrix, i, j-1)
        for pat in pats:
            pat[0] = pat[0] + s1
            pat[1] = pat[1] + s2
            pat[2] = pat[2] + char
            paths.append(pat)
    if "L" in matrix[i][j][3]:
        s1 = seq1[i-1]
        s2 = "-"
        char = " "
        pats = bactrackAligne(seq1, seq2, matrix,  scoreMatrix, i-1, j)
        for pat in pats:
            pat[0] = pat[0] + s1
            pat[1] = pat[1] + s2
            pat[2] = pat[2] + char
            paths.append(pat)
    return paths

def findMax(matrix):
    """
    Find the position of the bigger element in a matrix.
    return that position
    """
    maximum = 0
    pos = (0, 0)
    for i in range(len(matrix)-1, 0, -1):
        for j in range(len(matrix[i])-1, 0, -1):
            if matrix[i][j][0] > maximum:
                maximum = matrix[i][j][0]
                pos = (i, j)
    return pos

def localAlignement(seq1, seq2, openingPenality, extendingPenality, score, scoreName):
    """
    Make a local alignment between two sequences (seq1 ans seq2)
    """
    matrix = createMatrix(seq1, seq2, openingPenality, extendingPenality, score, True)
    print("Local alignment between:")
    print(seq1)
    print(seq2)
    print("Using penalty {}/{}".format(openingPenality, extendingPenality))
    print("And substitution matrix: {}\n".format(scoreName.name))
    pos = findMax(matrix)
    for i in bactrackAligne(seq1, seq2, matrix, score, pos[0], pos[1]):
        s1, s2, chars = i
        localScore = matrix[pos[0]][pos[1]][0]
        ident = chars.count(":") *100 / len(chars)
        sim = chars.count(".") *100 / len(chars) + ident
        print("-"*len(chars))
        print("Score: {}, Identity: {:.2f}%, similarity: {:.2f}%".format(localScore, ident, sim))
        print(s1)
        print(chars)
        print(s2)



def globalAlignement(seq1, seq2, openingPenality, extendingPenality, score, scoreName):
    """
    Make a global alignment between two sequences (seq1 ans seq2)
    """
    matrix = createMatrix(seq1, seq2, openingPenality, extendingPenality, score)
    globalScore = matrix[-1][-1][0]
    print("Global alignment between:")
    print(seq1)
    print(seq2)
    print("Using penalty {}/{}".format(openingPenality, extendingPenality))
    print("And substitution matrix: {}\n".format(scoreName.name))
    for i in bactrackAligne(seq1, seq2, matrix, score, len(seq1), len(seq2)):
        ident = i[2].count(":") *100 / len(i[2])
        sim = i[2].count(".") *100 / len(i[2]) + ident
        print("-"*len(i[2]))
        print("Score: {}, Identity: {:.2f}%, similarity: {:.2f}%".format(globalScore, ident, sim))
        print(i[0])
        print(i[2])
        print(i[1])


def main(args):
    seq = list(Sequence.fromFile(args.sequenceFilename))
    score = Score.fromFile(args.substitutionMatrix)
    openingPenality = args.openingPenality
    extendingPenality = args.extendingPenality
    local = args.local

    if not args.nseq1:
        for i in range(len(seq)):
            seq1 = seq[i]
            for j in range(i + 1, len(seq)):
                seq2 = seq[j]
                if local:
                    localAlignement(seq1, seq2, openingPenality, extendingPenality, score, args.substitutionMatrix)
                else:
                    globalAlignement(seq1, seq2, openingPenality, extendingPenality, score, args.substitutionMatrix)
                print("\n" + "#"*80 + "\n")
    else:
        if not (1 <= args.nseq1 <= len(seq) and 1 <= args.nseq2 <= len(seq)):
            parser.error("-1 and -2 option must be in range of the number of sequences.")
        seq1 = seq[args.nseq1 - 1]
        seq2 = seq[args.nseq2 - 1]
        if local:
            localAlignement(seq1, seq2, openingPenality, extendingPenality, score, args.substitutionMatrix)
        else:
            globalAlignement(seq1, seq2, openingPenality, extendingPenality, score, args.substitutionMatrix)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Make the global (or local) alignment between multiple sequences")
    parser.add_argument("sequenceFilename", type=argparse.FileType('r'), help="The file containing all the sequences to aligne.")
    parser.add_argument("substitutionMatrix", type=argparse.FileType('r'), help="The file containing the substitution matrix to use.")
    parser.add_argument("openingPenality", type=int, help="The penalty of opening a new gap.")
    parser.add_argument("extendingPenality", type=int, help="The penalty of extending a gap.")
    parser.add_argument("-l", "--local", action="store_true", help="Does a local alignment instead of a global.")
    parser.add_argument("-1", "--nseq1", type=int, help="The number of the first sequence to parse in file (if None, all sequences).")
    parser.add_argument("-2", "--nseq2", type=int, help="The number of the second sequence to align in file (required if -1 is used).")

    args = parser.parse_args()
    if args.nseq1 and not args.nseq2:
        parser.error("-2 option required if -1 option used.")
    elif args.nseq2 and not args.nseq1:
        parser.error("C'ant use -2 option if -1 is not used.")
    main(args)
