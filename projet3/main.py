
from math import log, sqrt
from sequence import Sequence
import argparse


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
    consensus = Sequence()
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

def printPSSM(PSSM):
    print("\n===== PSSM =====\n")
    key = PSSM[0].keys()
    for k in key:
        print("{:^7}".format(k), end=" ")
    print()
    for col in PSSM:
        for m in key:
            print("{: 7.3f}".format(col[m]), end=" ")
        print()


def main(args):
    sequenceFilename = args.sequenceFilename
    probaFilename = args.probaFilename
    occurences, nbSeq = getOccurences(sequenceFilename)

    proba = getProba(probaFilename)
    beta = sqrt(nbSeq)
    if args.beta:
        beta = args.beta

    PSSM, consensus = getPSSM(occurences, proba, beta, nbSeq)
    print("\n===== consensus =====\n")
    print(consensus)
    printPSSM(PSSM)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("sequenceFilename", help="The file containing all the sequences to parse.")
    parser.add_argument("probaFilename", help="The file containing the probability for all amino acid.")
    parser.add_argument("-b", "--beta", help="Set the value to use for beta (default = sqrt(nombre de sequences).", type=float)

    args = parser.parse_args()
    main(args)
