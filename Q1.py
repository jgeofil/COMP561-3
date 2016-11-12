import numpy as np

f =  open('test.fasta')
seq = ''.join(f.read().splitlines()[1:])


def hydroViterbi (Y):
    T = len(Y)
    S = ['O','X','M']
    K = len(S)
    B = {
        'O': {}
    }

    D1 = np.zeros((K,T))
    D2 = np.zeros((K,T))


hydroViterbi(seq)
