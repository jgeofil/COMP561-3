import numpy as np
import sys, getopt
import csv, math

def main(argv):
    helpMsg = 'H3Q1.py -i <inputfile>'
    inputfile = 'hw3_proteins.fa'

    try:
        opts, args = getopt.getopt(argv,"hi:")
    except getopt.GetoptError:
        print(helpMsg)
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print(helpMsg)
            sys.exit()
        elif opt == "-i":
            inputfile = arg

    print 'Opening file ' + inputfile +' (specify -i to use other file)...'

    f =  open(inputfile)
    lines = f.read().split('>')[1:]

    seqs = [''.join(l.splitlines()[1:]) for l in lines]
    names = [l.splitlines()[0] for l in lines]

    def hydroViterbi (Y):
        T = len(Y)
        S = ['I','O','M']
        K = len(S)
        aO = ['A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W']
        aI = ['R', 'H', 'K', 'D', 'E', 'S', 'T', 'N', 'Q', 'C', 'G', 'P']

        def Iprob (a):
            return (1/15.0 if a in aI else 1/40.0)
        def Oprob (a):
            return (3/40.0 if a in aO else 1/30.0)
        def Mprob (a):
            return 1/20.0

        A = np.matrix([
            [7/8.0, 3/80.0, 7/80.0],
            [1/25.0, 4/5.0, 4/25.0],
            [1/14.0, 1/14.0, 6/7.0]
        ])
        B = [Iprob, Oprob, Mprob]
        Init = [1/20.0]*len(S)

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
                        values.append(D1[k, i-1] + math.log(A[k,j]))
                    D1[j,i] = np.max(values) + math.log(B[j](a))
                    D2[j,i] = np.argmax(values)

        idx = np.argmax(D1[:,-1])
        path = S[idx]

        for i in range(T-1,0,-1):
            idx = D2[idx,i]
            path += S[idx]

        return path[::-1]

    print 'Running Viterbi algorithm...'
    cons = [hydroViterbi(s) for s in seqs]

    print 'Viterbi done.'
    print 'Wrtting output to H3Q1_out.fa...'
    f =  open('H3Q1_out.fa', 'w')
    for i,s in enumerate(cons):
        f.writelines(s+'\n')


    Xcount = [[len(ss) for ss in s.replace('I','M').split('M')] for s in cons]
    Xcount = [np.max(c) for c in Xcount]
    print '--------------------------------------------------------------------'
    print 'c.1 Which protein contains the longest hydrophobic region?'
    print 'Longest hydrophobic region ( Len = ', np.max(Xcount), ')'
    pos = np.argmax(Xcount)
    print names[pos]

    print '--------------------------------------------------------------------'
    print 'c.2 Which protein contains the largest fraction of amino acids annotated as belonging to a Mixed region?'
    Mcount = np.array([s.count('M')/float(len(s)) for s in cons])
    mlen = np.array([len(s) for s in cons])
    maxim = np.max(Mcount)
    res = np.array(names)[Mcount >= maxim]
    lens = mlen[Mcount >= maxim]
    relen = np.argmax(lens)
    print 'A total of ' + str(len(res)) + ' proteins with fraction = ' + str(maxim)
    print 'Of which the longest is:'
    print res[relen]

    print '--------------------------------------------------------------------'
    print 'c.4 What are the observed overall length distributions of Hydrophobic, Hydrophilic, and Mixed regions? '
    [item for sublist in l for item in sublist]
    Icount = np.array([len(ss) for s in cons for ss in s.replace('O','M').split('M') ])
    Ocount = np.array([len(ss) for s in cons for ss in s.replace('I','M').split('M') ])
    Mcount = np.array([len(ss) for s in cons for ss in s.replace('O','I').split('I') ])

    Icount = Icount[Icount != 0]
    Ocount = Ocount[Ocount != 0]
    Mcount = Mcount[Mcount != 0]

    '''

    with open('O.csv', 'wb') as csvfile:
        csvwrwite = csv.writer(csvfile, delimiter='\t',
                                quotechar='|', quoting=csv.QUOTE_MINIMAL)
        for i in Ocount:
            csvwrwite.writerow([i])
    with open('I.csv', 'wb') as csvfile:
        csvwrwite = csv.writer(csvfile, delimiter='\t',
                                quotechar='|', quoting=csv.QUOTE_MINIMAL)
        for i in Icount:
            csvwrwite.writerow([i])
    with open('M.csv', 'wb') as csvfile:
        csvwrwite = csv.writer(csvfile, delimiter='\t',
                                quotechar='|', quoting=csv.QUOTE_MINIMAL)
        for i in Mcount:
            csvwrwite.writerow([i])

    '''

    print 'Length counts written to: distributions.csv'
    print 'Avg(O) :' + str(np.mean(Ocount))
    print 'Avg(I) :' + str(np.mean(Icount))
    print 'Avg(M) :' + str(np.mean(Mcount))

    maxlen = max(np.amax(Icount)+1, np.amax(Ocount)+1, np.amax(Mcount)+1)
    Icount = np.bincount(Icount, minlength=maxlen)
    Ocount = np.bincount(Ocount, minlength=maxlen)
    Mcount = np.bincount(Mcount, minlength=maxlen)


    with open('distributions.csv', 'wb') as csvfile:
        csvwrwite = csv.writer(csvfile, delimiter='\t',
                                quotechar='|', quoting=csv.QUOTE_MINIMAL)
        csvwrwite.writerow(['Length','Hydrophobic', 'Hydrophilic', 'Mixed'])
        for i in range(0, maxlen):
            csvwrwite.writerow([i+1,Ocount[i], Icount[i], Mcount[i]])


    print '--------------------------------------------------------------------'
    print 'c.4 What are the observed amino acid frequencies in Hydrophobic, Hydrophilic, and Mixed regions?'
    aminos = ['A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W', 'R', 'H', 'K', 'D', 'E', 'S', 'T', 'N', 'Q', 'C', 'G', 'P']
    I,O,M = {},{},{}
    for am in aminos:
        I[am],O[am],M[am] = 0,0,0
    Itot = 0
    Otot = 0
    Mtot = 0
    for r,s in zip(cons, seqs):
        for x,a in zip(r,s):
            if x == 'O':
                Otot += 1
                O[a] += 1
            if x == 'I':
                Itot += 1
                I[a] += 1
            if x == 'M':
                Mtot += 1
                M[a] += 1
    for a in I:

        I[a] = I[a]/float(Itot)
    for a in O:
        O[a] = O[a]/float(Otot)
    for a in M:
        M[a] = M[a]/float(Mtot)

    with open('frequencies.csv', 'wb') as csvfile:
        print 'Amino acids frequencies written to: frequencies.csv'
        print ''
        csvwrwite = csv.writer(csvfile, delimiter='\t',
                                quotechar='|', quoting=csv.QUOTE_MINIMAL)
        csvwrwite.writerow(['Amino','Hydrophobic', 'Hydrophilic', 'Mixed'])
        print '\tHydrophobic\tHydrophilic\tMixed'
        for a in aminos:
            print a + ': ' + '\t'+ str(O[a]) + '\t'+ str(I[a]) + '\t'+ str(M[a])
            csvwrwite.writerow([a,O[a], I[a], M[a]])




if __name__ == "__main__":
   main(sys.argv[1:])
