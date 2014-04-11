from copy import deepcopy
from math import log

def fromFile(fileName):
    with open(fileName, 'r') as f:
        sequences = []
        for line in f:
            sequences.append(line.strip())
    return sequences

def similar(seq1, seq2, percentage):
    same = 0
    for i in range(min(len(seq1), len(seq2))):
        if(seq1[i] == seq2[i]):
            same += 1
    same = same / max(len(seq1), len(seq2))
    return same >= percentage

def findGroups(sequences, percentage):
    groups = []
    for sequence in sequences:
        if len(groups) == 0:
            groups.append([sequence])
        else:
            found = False
            for group in groups:
                for seq in group:
                    if similar(sequence, seq, percentage):
                        group.append(sequence)
                        found = True
                        break
                if found:
                    break
            else:
                groups.append([sequence])
    return groups

def calcGroup(i, groups):
    groupNb = 0
    pos = len(groups[groupNb])
    for group in groups:
        if i >= pos:
            groupNb += 1
            pos += len(groups[groupNb])
    return groupNb

def weightedGroup(i, groups):
    groupNb = 0
    pos = len(groups[groupNb])
    for group in groups:
        if i >= pos:
            groupNb += 1
            pos += len(groups[groupNb])
    return 1 / len(groups[groupNb])



def frequence(protein1, protein2, groups):
    summ = 0
    for i in range(len(groups[0][0])):
        col = map(lambda x: x[i], groups)
        for (i, prot) in enumerate(col):
            if prot == protein1:
                start = 0
                if protein1 == protein2:
                    start = i
                for j in range(start, len(col)):
                    if j != i and col[j] == protein2 and\
                       calcGroup(i, groups) != calcGroup(j, groups):
                        summ += weightedGroup(i, groups) *\
                                weightedGroup(j, groups)
    return summ

def weightedMatrix(groups, proteinList):
    matrix = []
    for i in range(len(proteinList)):
        matrix.append([0] * len(proteinList))
        for j in range(len(proteinList)):
            matrix[i][j] = frequence(i, j, groups)
    return matrix

def matrixAverage(matrixes):
    score = [[0 for i in range (20)] for j in range(20)]
    for i in range(len(matrixes[0])):
        for j in range(len(matrixes[0][0])):
            for n in range(len(matrixes)):

                score[i][j] += matrixes[n][i][j]
            score[i][j] = score[i][j] / len(matrixes)
    return score

def sumUpPartMatrix(matrix):
    summ = 0
    for i in range(len(matrix)):
        for j in range(len(matrix[i])):
            summ += matrix[i][j]
    return summ

def occurenceMatrix(matrix):
    sumUpPart = sumUpPartMatrix(matrix)
    if(sumUpPart !=0):
        for i in range(len(matrix)):
            for j in range(len(matrix[i])):
                matrix[i][j] = matrix[i][j] / sumUpPart
    return matrix

def residueFrequence(matrix, i):
    summ = 0
    for j in range(len(matrix)):
        if i != j:
            summ += matrix[j][i]
    summ /= 2
    summ += matrix[i][i]
    return summ


def alignementFrequence(matrix, i, j):
    if i == j:
        freq = residueFrequence(matrix, j) ** 2
    else:
        freq = 2 * residueFrequence(matrix, i) * residueFrequence(matrix, j)
    return freq

def logChance(matrix, i, j):
    try:
        if matrix[i][j] / alignementFrequence(matrix, i, j) == 0:
            res = 0
        else:
            res = 2 * log(matrix[i][j] / alignementFrequence(matrix, i, j), 2)
    except(ZeroDivisionError):
        res = 0
    return res

def calcEndMatrix(matrix):
    copyMatrix = deepcopy(matrix)
    for i in range(len(matrix)):
        for j in range(len(matrix[i])):
            matrix[i][j] = logChance(copyMatrix, i, j)
    return matrix

def main():
    matrixes = []
    fileList = ["PR00109_0.txt", "PR00109_1.txt", "PR00109_2.txt", "PR00109_3.txt", "PR00109_4.txt"]
    percentage = 0
    proteinList = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    for fileName in fileList:
        sequences = fromFile(fileName)
        groups = findGroups(sequences, percentage)
        matrixes.append(weightedMatrix(groups, proteinList))
    matrix = matrixAverage(matrixes)
    matrix = occurenceMatrix(matrix)
    matrix = calcEndMatrix(matrix)
    print(matrix)

if __name__=='__main__':
    main()

