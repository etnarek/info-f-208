#! /bin/python3

import sys
from sequence import Sequence
from score import Score

def createMatrix(seq1, seq2, I, E, score):
    matrix = [[[0, 0, 0, ""] for j in seq2] for i in seq1]
    matrix[0][0] = [0, -float("inf"), -float("inf"), ""]
    for i in range(1, len(seq1)):
        matrix[i][0] = [-float("inf"), -I - (i-1)*E, -float("inf"), ""]
    for i in range(1, len(seq2)):
        matrix[0][i] = [-float("inf"), -float("inf"), -I - (i-1)*E, ""]

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
            if S >= V or S >= W:
                matrix[i][j][3] = "D"
            elif V >= W:
                matrix[i][j][3] = "L"
            else:
                matrix[i][j][3] = "U"

            S = max(S, V, W)
            matrix[i][j][0] = S
    return matrix

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
            print(matrix)
            aligne(seq1, seq2, matrix, score, I, E)
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
