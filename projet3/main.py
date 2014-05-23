
from math import log, sqrt
from sequence import Sequence


def getOccurences(filename):
    occurences = []
    nbSeq = 0
    for sequence in Sequence.fromFile(filename):
        nbSeq += 1
        for i, aa in enumerate(sequence):
            if i < len(occurences):
                if aa != "-":
                    occurences[i][aa] = occurences[i].get(aa, 0) + 1
            else:
                if aa != "-":
                    occurences.append({aa: 1})
                else:
                    occurences.append({})
    return occurences, nbSeq


def getProba(filename):
    proba = {}
    with open(filename) as f:
        for l in f:
            if l[0] != "#":
                aa, prob = l.split(' ')
                proba[aa] = float(prob) / 100
    return proba


def getPSSM(occurences, proba, beta, nbSeq):
    PSSM = []
    consensus = []
    for i in range(len(occurences)):
        PSSM.append({})
        muaMax = -float("inf")
        aaMax = ""
        for aa in proba.keys():
            qua = (occurences[i].get(aa, 0) + beta * proba[aa])
            qua /= (nbSeq + beta)
            if qua == 0:
                mua = -float("inf")
            elif proba[aa] == 0:
                mua = float("inf")
            else:
                mua = log(qua / proba[aa], 10)
            PSSM[i][aa] = mua
            if mua > muaMax:
                muaMax = mua
                aaMax = aa
        consensus.append(aaMax)
    return PSSM, consensus


def main(sequenceFilename, probaFilename):
    occurences, nbSeq = getOccurences(sequenceFilename)
    print(occurences)
    print(nbSeq)

    proba = getProba(probaFilename)
    print(proba)
    beta = sqrt(nbSeq)
    beta = 1
    print(beta)

    PSSM, consensus = getPSSM(occurences, proba, beta, nbSeq)
    print(consensus)

    moy = 0
    maxm = -float("inf")
    i = 0
    key = PSSM[0].keys()
    for k in key:
        print(k, end=" ")
    print()
    for col in PSSM:
        for m in key:
            print("{:.4f}".format(col[m]), end=" ")
            if col[m] != float("inf") and col[m] != -float("inf"):
                moy += col[m]
                i += 1
                if col[m] > maxm:
                    maxm = col[m]
        print()

    print("moyenne: ", moy/i)
    print("max :", maxm)


if __name__ == "__main__":
    main("msaresults-MUSCLE.fasta", "probafile.txt")
