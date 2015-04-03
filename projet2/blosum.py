#!/usr/bin/python3
import sys
import argparse
from score import Score
from math import log

def blocksFromFiles(files):
    """load and return blocks from the given files
    """
    blocks = []
    for f in files:
        sequences = f.readlines()
        sequences = list(map(lambda x: x.strip(), sequences))
        blocks.append(sequences)
    return blocks

def haveEnoughtIdentity(seq1, seq2, percentage):
    """Calculate identity between two sequences
    return True if identity between seq1 and seq2 > percentage%
    """
    identical = 0
    different = 0
    for aa1, aa2 in zip(seq1, seq2):
        if aa1 == aa2:
            identical += 1
            if identical / float(len(seq1)) >= percentage / 100.0:
                return True
        else:
            different += 1
            if different / float(len(seq1)) > 1 - (percentage / 100.0):
                return False

    return identical / float(len(seq1)) >= percentage / 100.0

def canPlace(sequence, group, percentage):
    """ Verify if a sequence can be placed in a given group with a percentage
    """
    for seq in group:
        if haveEnoughtIdentity(seq, sequence, percentage):
            return True
    return False

def placeInGroups(sequences, percentage):
    """
    Make groups with all the sequences
    """
    groups = []
    for seq in sequences:
        for group in groups:
            if canPlace(seq, group, percentage):
                group.append(seq)
                break
        else:
            groups.append([seq])
    return groups

def ponderateFrequencyMatrix(groups):
    """
    Return the weighted matrix for the groups
    """
    matrix = Score()
    length = len(groups[0][0])
    for i in range(length):
        for iGroup, group in enumerate(groups):
            for seq in group:
                aa1 = seq[i]
                for iOtherGroup, otherGroup in enumerate(groups):
                    if iGroup != iOtherGroup:
                        for otherSeq in otherGroup:
                            aa2 = otherSeq[i]
                            freq = (1 / len(group)) * (1 / len(otherGroup))
                            if aa1 == aa2:
                                freq = freq / 2
                            matrix[aa1, aa2] = matrix.get((aa1, aa2), 0) + freq
    return matrix

def calculateFreqSum(matrix):
    """
    Calculate the sum of the frequencies of a matrix.
    """
    diag, tot = 0, 0
    for key, val in matrix.items():
        if key[0] == key[1]:
            diag += val
        else:
            tot += val
    tot = tot / 2 + diag
    return tot

def occurenceProbability(matrix):
    """
    Calculate the occurrence weighted matrix of a frequencies matrix
    """
    ponderate = Score()
    freqSum = calculateFreqSum(matrix)

    for key, val in matrix.items():
        ponderate[key] = val / freqSum
    return ponderate

def propability(matrix, aa):
    """
    give the probability of a replacement for an aa
    """
    prob = 0
    for key, val in matrix.items():
        if key[0] == aa and key[1] != aa:
            prob += val / 2
    return prob + matrix.get((aa, aa), 0)

def expectedFrequency(matrix):
    """
    give the evolution matrix
    """
    expected = Score()

    for key in matrix.keys():
        if key[0] == key[1]:
            expected[key] = propability(matrix, key[0]) ** 2
        else:
            expected[key] = 2 * propability(matrix, key[0]) * propability(matrix, key[1])
    return expected

def loggOddRatio(ponderate, expected, blosum):
    """
    return the log odd ratio matrix.
    """
    for key in ponderate.keys():
        blosum[key] = blosum.get(key, 0) + 2 * log(ponderate[key] / expected[key], 2)
    return blosum

def average(blosum, nbBlocks):
    """
    Compute the average of the matrix by the number of blocks.
    """
    for key in blosum.keys():
        blosum[key] = round(blosum[key] / nbBlocks)
    return blosum

def main(args):
    percentage = args.identity
    blosum = Score()
    blocks = blocksFromFiles(args.blocks)

    for sequences in blocks:
        groups = placeInGroups(sequences, percentage)
        frequencies = ponderateFrequencyMatrix(groups)
        ponderate = occurenceProbability(frequencies)
        expected = expectedFrequency(ponderate)
        blosum = loggOddRatio(ponderate, expected, blosum)
    blosum = average(blosum, len(blocks))
    print(blosum)
    if args.output:
        args.output.write(str(blosum))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create blosum matrix from blocs")
    parser.add_argument("identity", type=int, help="The percentage of identity that 2 sequences must have.")
    parser.add_argument('blocks', type=argparse.FileType('r'), nargs='+', help="All the files each containing a block to use to create the blosum matrix.")
    parser.add_argument('-o', "--output", type=argparse.FileType('w'), help="The file to write the generated matrix.")

    args = parser.parse_args()
    main(args)

