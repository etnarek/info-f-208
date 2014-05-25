import argparse


def main(args):
    pass


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate the profile (PSSM matrix) with the sequences in the sequenceFile and the probability in the probaFile. The output can be set in the outFile.")
    parser.add_argument("infoFile", help="The file containing all filename to parse.")
    parser.add_argument("probaFilename", help="The file containing the probability for all amino acid.")
    parser.add_argument("-b", "--beta", help="Set the value to use for beta (default = sqrt(nombre de sequences).", type=float)
    parser.add_argument("-o", "--outFile", help="The file to save the PSSM matrix.")

    args = parser.parse_args()
    main(args)
