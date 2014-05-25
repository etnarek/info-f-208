def parseDssp(filename, tmpFastaFile, folderName):
    with open(filename) as f:
        with open(tmpFastaFile, 'w') as tmpf:
            for l in f:
                treat(l.split()[0], tmpf, folderName)


def treat(line, tmpf, folderName):
    chaine = line[4]
    if folderName[-1] != "/":
        folderName += "/"
    filename = folderName + line[:4] + ".dssp"
    tmpf.write("> " + filename)

    with open(filename) as f:
        ll = f.readlines()
        tmpf.write(" | " + getName(ll[3]))
        tmpf.write(" | " + getOrganism(ll[4]) + "\n")
        seqs = getSeqs(ll[28:], chaine)
        tmpf.write(seqs[0] + "\n" + seqs[1] + "\n")


def getName(line):
    line = line.split()[3:]
    return " ".join(line)[:-1]


def getOrganism(line):
    line = line.split()[3:]
    return " ".join(line)[:-1]


def getSeqs(lines, chaine):
    resseq1 = ""
    resseq2 = ""
    for l in lines:
        if l[11] == chaine:
            resseq1 += l[13].upper()
            if l[16] in "HGIETCSB":
                resseq2 += family(l[16])
            else:
                resseq2 += "C"
    return (resseq1, resseq2)


def family(char):
    if char == 'H' or char == 'G' or char == 'I':
        return 'H'
    elif char == 'C' or char == 'S' or char == 'B':
        return 'C'
    else:
        return char
