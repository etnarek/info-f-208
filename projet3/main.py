
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
    exclude = []
    with open(filename) as f:
        for l in f:
            if l[0] != "#":
                aa, prob = l.split(' ')
                if float(prob) != 0:
                    proba[aa] = float(prob) / 100
                else:
                    exclude.append(aa)
    return proba, exclude


def getPSSM(occurences, proba, exclude, beta, nbSeq):
    PSSM = []
    consensus = Sequence()
    for i in range(len(occurences)):
        PSSM.append({})
        muaMax = -float("inf")
        aaMax = "-"
        for aa in proba.keys():
            if aa not in exclude:
                qua = (occurences[i].get(aa, 0) + beta * proba[aa])
                qua /= (nbSeq + beta)
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
    print("     ", end=" ")
    for k in key:
        print("{:^7}".format(k), end=" ")
    print()
    for i, col in enumerate(PSSM):
        print("{:5}".format(i), end=" ")
        for m in key:
            print("{: 7.3f}".format(col[m]), end=" ")
        print()


def savePSSM(PSSM, outFilename):
    with open(outFilename, 'w') as f:
        key = PSSM[0].keys()
        f.write("   ")
        for k in range(len(PSSM)):
            f.write("{:^7} ".format(k))
        f.write("\n")
        for m in key:
            f.write("{:2} ".format(m))
            for i, col in enumerate(PSSM):
                f.write("{: 7.3f} ".format(col[m]))
            f.write("\n")


def main(args):
    sequenceFilename = args.sequenceFilename
    probaFilename = args.probaFilename
    occurences, nbSeq = getOccurences(sequenceFilename)

    proba, exclude = getProba(probaFilename)
    beta = sqrt(nbSeq)
    if args.beta:
        beta = args.beta

    PSSM, consensus = getPSSM(occurences, proba, exclude, beta, nbSeq)
    print("\n===== consensus =====\n")
    print(consensus)
    printPSSM(PSSM)

    if args.outFile:
        savePSSM(PSSM, args.outFile)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate the profile (PSSM matrix) with the sequences in the sequenceFile and the probability in the probaFile. The output can be set in the outFile.")
    parser.add_argument("sequenceFilename", help="The file containing all the sequences to parse.")
    parser.add_argument("probaFilename", help="The file containing the probability for all amino acid.")
    parser.add_argument("-b", "--beta", help="Set the value to use for beta (default = sqrt(nombre de sequences).", type=float)
    parser.add_argument("-o", "--outFile", help="The file to save the PSSM matrix.")

    args = parser.parse_args()
    main(args)
