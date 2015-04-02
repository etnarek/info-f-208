import sys

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
    for seq in group:
        if haveEnoughtIdentity(seq, sequence, percentage):
            return True
    return False

def placeInGroups(sequences, percentage):
    groups = []
    for seq in sequences:
        for group in groups:
            if canPlace(seq, group, percentage):
                group.append(seq)
                break
        else:
            groups.append([seq])
    return groups

def main():
    sequences = ["ATCKQ", "SSCRN", "ATCRN", "TECRQ", "ASCKN", "SECEN", "SDCEQ"]
    percentage = 50

    groups = placeInGroups(sequences, percentage)
    print(groups)

if __name__ == "__main__":
    main()

