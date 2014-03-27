def createF(seq1, seq2, penality, score):
    F = [[0 for j in seq2] for i in seq1]
    for i in range(len(seq1)):
        F[i][0] = i * penality
    for i in range(len(seq2)):
        F[0][i] = i * penality

    for i in range(1, len(seq1)):
        for j in range(1, len(seq2)):
            choice1 = F[i-1][j-1] + score[seq1[i], seq2[j]]
            choice2 = F[i-1][j] + penality
            choice3 = F[i][j - 1] + penality
            F[i][j] = max(choice1, choice2, choice3)

    return F
