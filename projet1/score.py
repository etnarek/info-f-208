class Score(dict):
    def __int__(self, *arg, **kw):
        dict.__init__(*arg, **kw)

    @staticmethod
    def fromFile(filename):
        matrice = Score()
        with open(filename) as f:
            for line in f:
                if line[0] == ' ':
                    up = line.split()
                elif line[0] != '#':
                    head = line[0]
                    data = line[1:].split()
                    for i in range(len(data)):
                        matrice[up[i], head] = data[i]
        return matrice

