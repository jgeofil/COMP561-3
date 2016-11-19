import numpy as np

f =  open('hw3_proteins.fa')
lines = f.read().split('>')[1:]

seqs = [''.join(l.splitlines()[1:]) for l in lines]
names = [l.splitlines()[0][:10] for l in lines]

def hydroViterbi (Y):
    T = len(Y)
    S = ['I','O','M']
    K = len(S)
    aX = ['A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W']
    aO = ['R', 'H', 'K', 'D', 'E', 'S', 'T', 'N', 'Q', 'S', 'C', 'G', 'P']

    def Oprob (a):
        return (0.8/len(aO) if a in aO else 0.2/len(aX))
    def Xprob (a):
        return (0.6/len(aX) if a in aX else 0.4/len(aO))
    def Mprob (a):
        return 1/20.0

    A = np.matrix([
        [7/8.0, 3/80.0, 7/80.0],
        [1/25.0, 4/5.0, 4/25.0],
        [1/14.0, 1/14.0, 6/7.0]
    ])
    B = [Oprob, Xprob, Mprob]
    Init = [1.0/len(S)]*len(S)

    D1 = np.zeros((K,T))
    D2 = np.zeros((K,T),dtype=int)

    for i, s in enumerate(S):
        D1[i,0] = Init[i] * B[i](Y[0])
        D2[i,0] = 0

    for i, a in enumerate(Y):
        if i > 0:
            for j in range(K):
                values = []
                for k in range(K):
                    values.append(D1[k, i-1] * A[k,j])

                D1[j,i] = np.max(values) * B[j](a) * 10
                D2[j,i] = np.argmax(values)

    idx = np.argmax(D1[:,-1])
    path = S[idx]

    for i in range(T-1,0,-1):
        idx = D2[idx,i]
        path += S[idx]

    return path[::-1]


cons = [hydroViterbi(s) for s in seqs]

f =  open('out.fa', 'w')
for s in cons:
    f.writelines(s+'\n')

print cons

print '------------------'
print 'Longest hydrophobic region'
Xcount = [[len(ss) for ss in s.replace('O','M').split('M')] for s in cons]
Xcount = [np.max(c) for c in Xcount]
print 'Len = ', np.max(Xcount)
pos = np.argmax(Xcount)
print names[pos]

print '------------------'
print 'Proportion of mixed aminos'
Mcount = [s.count('M')/float(len(s)) for s in cons]
print Mcount
