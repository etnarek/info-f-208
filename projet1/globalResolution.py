#! /bin/python3

import sys
from sequence import Sequence
from score import Score

def createMatrix(seq1, seq2, openingPenality, extendingPenality, scoreMatrix, local = False):
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

def bactrackAligne(seq1, seq2, matrix, scoreMatrix, i=-1, j=-1):
    paths = []
    if "E" in matrix[i][j][3]:
        return [["", "", ""]]
    if "D" in matrix[i][j][3]:
        s1 = seq1[i]
        s2 = seq2[j]
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
        s2 = seq2[j]
        char = " "
        pats = bactrackAligne(seq1, seq2, matrix, scoreMatrix, i, j-1)
        for pat in pats:
            pat[0] = pat[0] + s1
            pat[1] = pat[1] + s2
            pat[2] = pat[2] + char
            paths.append(pat)
    if "L" in matrix[i][j][3]:
        s1 = seq1[i]
        s2 = "-"
        char = " "
        pats = bactrackAligne(seq1, seq2, matrix,  scoreMatrix, i-1, j)
        for pat in pats:
            pat[0] = pat[0] + s1
            pat[1] = pat[1] + s2
            pat[2] = pat[2] + char
            paths.append(pat)
    return paths


def aligne(seq1, seq2, F, S, I, E):
    align1 = ""
    align2 = ""
    alignChar = ""
    i = len(seq1) - 1
    j = len(seq2) - 1
    identity = 0
    globalScore = 0
    while i > 0 and j > 0:
        if(F[i][j][3] == "D"):
            align1 = seq1[i] + align1
            align2 = seq2[j] + align2
            if seq1[i] == seq2[j]:
                alignChar = ":" + alignChar
                identity += 1
            elif S[seq1[i], seq2[j]] >= 0:
                alignChar = "." + alignChar
            else:
                alignChar = " " + alignChar
            i -= 1
            j -= 1

        elif F[i][j][3] == "L":
            align1 = seq1[i] + align1
            align2 = "-" + align2
            alignChar = " " + alignChar
            i -= 1

        elif F[i][j][3] == "U":
            align1 = "-" + align1
            align2 = seq2[j] + align2
            alignChar = " " + alignChar
            j -= 1

    if i == 0 and j == 0:
        align1 = seq1[i] + align1
        align2 = seq2[j] + align2
        if seq1[i] == seq2[j]:
            alignChar = ":" + alignChar
            identity += 1
        elif S[seq1[i], seq2[j]] >= 0:
            alignChar = "." + alignChar
        else:
            alignChar = " " + alignChar
        i -= 1
        j -= 1

    while i >= 0:
        align1 = seq1[i] + align1
        align2 = "-" + align2
        i -= 1
    while j >= 0:
        align1 = "-" + align1
        align2 = seq2[j] + align2
        j -= 1
    globalScore = F[len(F)-1][len(F[0])-1][0]
    print(align1)
    print(alignChar)
    print(align2)
    print("% d'identité: {}".format(identity * 100/max(len(seq1), len(seq2))))
    print("score global: {}".format(globalScore))
    return align1, align2, alignChar, identity, globalScore


def main():
    fileName = sys.argv[1]
    seq = list(Sequence.fromFile(fileName))
    fileName = sys.argv[2]
    score = Score.fromFile(fileName)
    I = int(sys.argv[3])
    E = int(sys.argv[4])

    for i in range(1):
        seq1 = seq[i]
        for j in range(i + 1, 2):
            seq2 = seq[j]
            matrix = createMatrix(seq1, seq2, I, E, score)
            #print(matrix)
            #aligne(seq1, seq2, matrix, score, I, E)
            print(bactrackAligne(seq1, seq2, matrix, score))
            print()

if __name__ == '__main__':
    if len(sys.argv) == 5:
        main()
    else:
        print("Pour utiliser ce programme, il faut donner les information suivante:")
        print("Le fichier contenant les sequences")
        print("Le fichier contenant la matrice de substitution")
        print("La valeur d'initialisation d'un gap")
        print("La valeur pour continuer un gap")
        print("\nExample d'utilisation:")
        print("> python3 globalResolution.py data/PDZ-sequences.fasta data/blosum62.txt 14 4")
