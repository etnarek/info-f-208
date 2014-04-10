def createF(seq1, seq2, I, E, score):
    F = [[0 for j in seq2] for i in seq1]
    for i in range(len(seq1)):
        F[i][0] = i * I
    for i in range(len(seq2)):
        F[0][i] = i * I

    gap2 = 0
    gap3 = 0
    for i in range(1, len(seq1)):
        for j in range(1, len(seq2)):
            penality2 = getPenality(gap2, I, E)
            penality3 = getPenality(gap3, I, E)
            choice1 = F[i-1][j-1] + score[seq1[i], seq2[j]]
            choice2 = F[i-1][j] + penality2
            choice3 = F[i][j - 1] + penality3
            F[i][j] = max(choice1, choice2, choice3)
            if choice2 > choice1:
                penality2 += 1
                penality3 = 0
            elif choice3 > choice2:
                penality3 += 1
                penality2 = 0
            else:
                penality2, penality3 = 0, 0


    return F

def getPenality(nbGap, I, E):
    return -I-nbGap*E

def aligne(seq1, seq2, F):
    pass
