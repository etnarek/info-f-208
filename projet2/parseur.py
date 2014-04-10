import sys

def div(filename):
    blocks = []
    with open(filename, 'r') as f:
        if f.readline()[0] == '>':
            pass
        for line in f:
            if line[0] == '>':
                break
            blocks.append([line.strip()])
        i = 0
        for line in f:
            if line[0] == '>':
                i = 0
            else:
                blocks[i].append(line.strip())
                i += 1
    return blocks

def writeBlocks(blocks, filename):
    for (i, block) in enumerate(blocks):
        newFilename = filename.strip(".txt")
        newFilename += "_{}.txt".format(i)
        with open(newFilename, 'w') as f:
            for line in block:
                f.write("%s\n" % line)

def main():
    print("Hello")
    filename = sys.argv[1]
    blocks = div(filename)
    writeBlocks(blocks, filename)

if __name__=='__main__':
    main()
