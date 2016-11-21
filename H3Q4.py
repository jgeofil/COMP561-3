import math
import numpy as np
import sys, getopt
import scipy.stats as stats

K_LEN = 6

def openSequences (fileName):
    f =  open(fileName)
    lines = f.read().split('>')[1:]
    seqs = [''.join(l.splitlines()[1:]) for l in lines]
    return seqs

def recKmers (s, kmers, L):
    if len(s) >= L:
        kmers.append(s)
    else:
        for c in ['A', 'T', 'G', 'C', 'X', 'Y', 'Z']:
            recKmers(s + c, kmers, L)

def listKmers (L):
    kmers = []
    recKmers('', kmers, L)
    return kmers

class Tree(object):
    def __init__(self):
        self.children = []
        self.amino = None
        self.count = 0

    def addCount (self, kmer):
        if len(kmer) == 0:
            self.count += 1
        else:
            c = kmer[0]
            found = False
            for node in self.children:
                if node.amino == c:
                    found = True
                    node.addCount(kmer[1:])
            if not found:
                newNode = Tree()
                newNode.amino = c
                self.children.append(newNode)
                newNode.addCount(kmer[1:])

    def readSequences(self, seqs):
        for s in seqs:
            for i in range(len(s) - K_LEN):
                self.addCount(s[i:i+K_LEN])

    def getCount (self, kmer):
        if len(kmer) == 0:
            return self.count
        else:
            c = kmer[0]
            total = 0
            if c == 'X':
                c = ['A', 'G']
            elif c == 'Y':
                c = ['C', 'T']
            elif c == 'Z':
                c = ['A', 'G', 'C', 'T']
            else:
                c = [c]
            for node in self.children:
                if node.amino in c:
                    total += node.getCount(kmer[1:])
            return total

    def getCountsKmers (self, kmers):
        counts = []
        for k in kmers:
            counts.append(self.getCount(k))
        return np.array(counts)



PROBS = {'A': 0.237, 'T': 0.238, 'G': 0.261, 'C': 0.263, 'X': 0.237+0.261, 'Y': 0.238+0.263, 'Z': 1}

def getExpect (seqs, kmers):
    sumlen = sum([len(s) - K_LEN + 1 for s in seqs])
    expects = [reduce(lambda x, y: y*x, [PROBS[p] for p in k])*sumlen for k in kmers]
    return expects

def sumlen (seqs):
    sumlen = sum([len(s) - K_LEN + 1 for s in seqs])
    return float(sumlen)

def getZScore(count, notcount, Ecount, Enots, lenCount, lenNot):
    return [(Nw/float(Nn) - lenCount/float(lenNot)/math.sqrt(lenCount/float(lenNot))) for (Nw, Nn, Ew, En) in zip(count, notcount, Ecount, Enots)]

def main(argv):
    helpMsg = 'H3Q1.py -i <inputfile> -c <controlfile>'
    inputfile = 'GATA2_chr1.fa'
    controlfile = 'not_GATA2_chr1.fa'

    try:
        opts, args = getopt.getopt(argv,"hi:c:")
    except getopt.GetoptError:
        print(helpMsg)
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print(helpMsg)
            sys.exit()
        elif opt == "-i":
            inputfile = arg
        elif opt == "-c":
            controlfile = arg

    seqs = openSequences(inputfile)
    nots = openSequences(controlfile)

    print '--------------------------------------------------------------------'
    print 'Bound sequences: ' + str(len(seqs))
    print 'Unbound sequences: ' + str(len(nots))

    print (sumlen(nots)/sumlen(seqs))

    print 'Listing k-mers of length ' + str(K_LEN) + '...'
    kmers = listKmers(K_LEN)

    t = Tree()
    nt = Tree()

    print 'Counting k-mers in bound sequences...'
    t.readSequences(seqs)
    counts = t.getCountsKmers(kmers)

    print 'Counting k-mers in unbound sequences...'
    nt.readSequences(nots)
    ncounts = nt.getCountsKmers(kmers)

    print 'Calculating Z-score...'
    Eseqs = getExpect(seqs, kmers)
    Enots = getExpect(nots, kmers)
    #zPos = getZScore(counts, Eseqs)
    #zNeg = getZScore(ncounts, Enots)
    #norm = getNormZ(zPos, zNeg)

    norm = getZScore(counts, ncounts, Eseqs, Enots, sumlen(seqs), sumlen(nots))


    ma = np.argsort(norm)

    kmers = np.array(kmers)[ma]
    norm = np.array(norm)[ma]
    print '--------------------------------------------------------------------'
    print 'Motif\tZ-score\t\tX=A|G\tY=C|T\tZ=A|T|G|C'
    print ''
    for x in range(-1, -10, -1):
        print kmers[x] + '\t' + str(norm[x])
    for x in range(0, 10, 1):
        print kmers[x] + '\t' + str(norm[x])



if __name__ == "__main__":
   main(sys.argv[1:])
