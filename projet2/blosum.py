from score import Score
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
                if not found:
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


def colGroups(groups, i):
    col = []
    for group in groups:
        tmp = []
        for seq in group:
            tmp.append(seq[i])
        col.append(tmp)
    return col

def frequence(protein1, protein2, groups):
    summ = 0
    print("len :", len(groups[0][0]))
    for k in range(len(groups[0][0])):
        #print("k: ", k)
        col = colGroups(groups, k)
        for i in col:
            for j in col:
                if i != j:
                    summ += i.count(protein1)/len(i) *\
                            j.count(protein2)/len(j)
    return summ

def weightedMatrix(groups, proteinList):
    matrix = Score()
    for i in proteinList:
        for j in proteinList:
            matrix[i, j] = frequence(i, j, groups)
    return matrix

def matrixAverage(matrixes):
    score = Score()
    for matrix in matrixes:
        for i, j in matrix.keys():
            print(i, j)
            score[i, j] = matrix[i, j] + score.get((i, j), 0)
    for i, j in score.keys():
        score[i, j] /= len(matrixes)
    return score

def sumUpPartMatrix(matrix):
    summ = 0
    for i, j in matrix.keys():
        summ += matrix[i, j]
    return summ

def occurenceMatrix(matrix):
    sumUpPart = sumUpPartMatrix(matrix)
    if(sumUpPart !=0):
        for i, j in matrix.keys():
            matrix[i, j] = matrix[i, j] / sumUpPart
    return matrix

def residueFrequence(matrix, i):
    summ = 0
    keys = matrix.keys()
    keys_list = []
    keys_list = [x[0] for x in keys if x[0] not in keys_list]
    for j in keys_list:
        if i != j:
            summ += matrix[j, i]
    summ /= 2
    summ += matrix[i, i]
    return summ


def alignementFrequence(matrix, i, j):
    if i == j:
        freq = residueFrequence(matrix, j) ** 2
    else:
        freq = 2 * residueFrequence(matrix, i) * residueFrequence(matrix, j)
    return freq

def logChance(matrix, i, j):
    try:
        if matrix[i, j] / alignementFrequence(matrix, i, j) == 0:
            res = 0
        else:
            res = 2 * log(matrix[i, j] / alignementFrequence(matrix, i, j), 2)
    except(ZeroDivisionError):
        res = 0
    return res

def calcEndMatrix(matrix):
    copyMatrix = deepcopy(matrix)
    for i, j in matrix.keys():
        matrix[i, j] = logChance(copyMatrix, i, j)
    return matrix

def afficherMatrice(matrice):

    indiceLigne = 0

    while (indiceLigne != len(matrice)):
        indiceColonne = 0
        while (indiceColonne != len(matrice)):
            print(round(matrice[indiceLigne][indiceColonne]), end=' ' *
                  (4 - len(str(round(matrice[indiceLigne][indiceColonne])))))
            indiceColonne = indiceColonne + 1
        indiceLigne = indiceLigne + 1
        print()

def main():
    matrixes = []
    fileList = ["PR00109_0.txt", "PR00109_1.txt", "PR00109_2.txt", "PR00109_3.txt", "PR00109_4.txt"]
    percentage = 0.7
    proteinList = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    for fileName in fileList:
        sequences = fromFile(fileName)
        groups = findGroups(sequences, percentage)
        for group in groups:
            print(len(group))
        matrixes.append(weightedMatrix(groups, proteinList))
        print(0)
    print(1)
    matrix = matrixAverage(matrixes)
    print(2)
    matrix = occurenceMatrix(matrix)
    print(3)
    matrix = calcEndMatrix(matrix)
    print(4)
    afficherMatrice(matrix)

if __name__=='__main__':
    main()

