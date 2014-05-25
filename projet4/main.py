import argparse
from math import log,  sqrt
from parser import parseDssp

TMP_FASTA = "tmp.fasta"
FAMILIES = "HETC"
AA = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z"]


def main(args):
    filename = args.filename
    foldername = args.foldername

    print("----- Parse dssp files -----")
    parseDssp(filename, TMP_FASTA, foldername)
    print("----- Start GORIII algo -----")
    gor3(filename, foldername)


def gor3(testFilename, testFoldername):
    mAa, mSec, mVois, mIndiv, mPairs = initM()
    with open(TMP_FASTA) as tFasta:
        lines = tFasta.readlines()
        mSeq, mStruct = seqMatrixes(lines)
        print("      Calculate frequencies")
        fillFreq(mSeq, mStruct, mAa, mSec, mVois)
        print("      Calculate score matrix")
        score(mAa, mSec, mVois, mIndiv, mPairs)
        predict = [0 for i in range(len(mSeq))]
        Q3 = [0 for i in range(len(mSeq))]
        print("\n===== result =====")
        calcInfo(Q3, predict, mSeq, mStruct, mIndiv, mPairs)
        print("\n===== Total =====")
        Q3Tot(Q3)
        mccTot(predict, mSeq, mStruct)


def initM():
    fam = range(len(FAMILIES))
    aa = range(len(AA))
    mAa = [0 for i in fam]
    mSec = [[0 for i in aa] for j in fam]
    mVois = [[[0 for i in aa] for j in aa] for k in fam]
    mIndiv = [[0 for i in aa] for j in fam]
    mPairs = [[[0 for i in aa] for j in aa] for k in fam]
    return mAa, mSec, mVois, mIndiv, mPairs


def seqMatrixes(lines):
    mSeq = []
    mStruct = []
    for i, l in enumerate(lines):
        if i % 3 == 1:
            mSeq.append(l.strip())
        elif i % 3 == 2:
            mStruct.append(l.strip())
    return mSeq, mStruct


def fillFreq(mSeq, mStruct, mAa, mSec, mVois):
    for i in range(len(mSeq)):
        for j in range(len(mSeq[i])):
            aa = mSeq[i][j]
            struct = mStruct[i][j]
            s = FAMILIES.index(struct)
            r = AA.index(aa)
            mAa[s] += 1
            mSec[s][r] += 1
            for k in range(-8, 9):
                if k != 0:
                    indVois = j + k
                    if 0 <= indVois < len(mSeq[i]):
                        aav = mSeq[i][indVois]
                        rv = AA.index(aav)
                        mVois[s][rv][r] += 1


def score(mAa, mSec, mVois, mIndiv, mPairs):
    for i in range(len(mIndiv)):
        for j in range(len(mIndiv[0])):
            calcIndiv(i, j, mAa, mSec, mIndiv)
            for k in range(len(mPairs[0][0])):
                calcPairs(i, j, k, mSec, mVois, mPairs)


def calcIndiv(s, r, mAa, mSec, mIndiv):
    nOnly = 0
    nIndiv = 0
    for i in range(len(FAMILIES)):
        if i != s:
            nOnly += mAa[i]
            nIndiv += mSec[i][r]

    if mAa[s] == 0 or nIndiv == 0 or mSec[s][r] == 0:
        mIndiv[s][r] = -float("inf")
    else:
        mIndiv[s][r] = log(mSec[s][r]/nIndiv) + log(nOnly/mAa[s])


def calcPairs(s, r, rv, mSec, mVois, mPairs):
    nPairs = 0
    nIndiv = 0
    for i in range(len(FAMILIES)):
        if i != s:
            nPairs += mVois[i][rv][r]
            nIndiv += mSec[i][r]

    if mSec[s][r] == 0 or nPairs == 0 or mVois[s][rv][r] == 0:
        mPairs[s][rv][r] = -float("inf")
    else:
        mPairs[s][rv][r] = log(mVois[s][rv][r]/nPairs) + log(nIndiv/mSec[s][r])


def calcInfo(Q3, predict, mSeq, mStruct, mIndiv, mPairs):
    for i, seq in enumerate(mSeq):
        res = ""
        for j, aa in enumerate(seq):
            r = AA.index(aa)
            tmp = []
            iMax = 0
            for s in range(len(FAMILIES)):
                tmp.append(mIndiv[s][r])
                for k in range(-8, 9):
                    if k != 0:
                        indVois = j+k
                        if 0 <= indVois < len(seq):
                            aav = seq[indVois]
                            rv = AA.index(aav)
                            tmp[s] += mPairs[s][rv][r]
                if tmp[s] >= tmp[iMax]:
                    iMax = s
            res += FAMILIES[iMax]
        predict[i] = res
        show(predict, Q3, mSeq, mStruct, i)


def show(predict, Q3, mSeq, mStruct, i):
    print("\n----- Seq n°{} -----".format(i + 1))
    print("Prediction: ", predict[i])
    print("Real: ", mStruct[i])
    print("Verification of the quality of the estimations.")
    calcQ3(Q3, predict, mSeq, mStruct, i)
    calcMCC(predict, mSeq, mStruct, i)


def calcQ3(Q3, predict, mSeq, mStruct, i):
    corr = 0
    for j in range(len(mSeq[i])):
        if predict[i][j] == mStruct[i][j]:
            corr += 1
    res = corr / len(mSeq[i])
    Q3[i] = res
    print("Q3: ", res)


def calcMCC(predict, mSeq, mStruct, i):
    for s in range(len(FAMILIES)):
        struct = FAMILIES[s]
        tp, tn, fp, fn = 0, 0, 0, 0
        for j in range(len(mSeq[i])):
            if predict[i][j] == mStruct[i][j]:
                if predict[i][j] == struct:
                    tp += 1
                else:
                    tn += 1
            else:
                if mStruct[i][j] == struct:
                    fn += 1
                elif predict[i][j] == struct:
                    fp += 1
        den = (tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)
        if den != 0:
            mcc = (tp*tn - fp*fn)/sqrt(den)
        else:
            mcc = "Error"
        print("MCC{}: {}".format(FAMILIES[s], mcc))


def Q3Tot(Q3):
    length = len(Q3)
    tot = sum(Q3)
    print("Tot Q3: {}".format(tot/length))


def mccTot(predict, mSeq, mStruct):
    for s in range(len(FAMILIES)):
        struct = FAMILIES[s]
        tp, tn, fp, fn = 0, 0, 0, 0
        for i in range(len(mSeq)):
            for j in range(len(mSeq[i])):
                if predict[i][j] == mStruct[i][j]:
                    if predict[i][j] == struct:
                        tp += 1
                    else:
                        tn += 1
                else:
                    if mStruct[i][j] == struct:
                        fn += 1
                    elif predict[i][j] == struct:
                        fp += 1
        den = (tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)
        if den != 0:
            mcc = (tp*tn - fp*fn)/sqrt(den)
        else:
            mcc = "Error"
        print("Tot MCC{}: {}".format(FAMILIES[s], mcc))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Prediction of the secondary structure with the GORIII algorithm.")
    parser.add_argument("filename", help="The file containing the sequences to parse.")
    parser.add_argument("foldername", help="The folder containing the data.")

    args = parser.parse_args()
    main(args)
