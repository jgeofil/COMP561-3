import math
import numpy as np

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
        return counts



PROBS = {'A': 0.25, 'T': 0.25, 'G': 0.25, 'C': 0.25, 'X': 0.5, 'Y': 0.5, 'Z': 1}

def getExpect (seqs, kmers):
    sumlen = sum([len(s) - K_LEN + 1 for s in seqs])
    expects = [reduce(lambda x, y: y*x, [PROBS[p] for p in k])*sumlen for k in kmers]
    return expects

def getZScore(counts, expects):
    return [(Nw - Ew)/math.sqrt(Ew) for (Nw,Ew) in zip(counts, expects)]

def getNormZ(zPos, zNeg):
    return [p/n for (p,n) in zip(zPos, zNeg)]

seqs = openSequences('GATA2_chr1.fa')
nots = openSequences('not_GATA2_chr1.fa')

kmers = listKmers(6)
print len(kmers)

t = Tree()
nt = Tree()

t.readSequences(seqs)
nt.readSequences(nots)

Eseqs = getExpect(seqs, kmers)
Enots = getExpect(nots, kmers)

counts = t.getCountsKmers(kmers)
ncounts = nt.getCountsKmers(kmers)

zPos = getZScore(counts, Eseqs)
zNeg = getZScore(ncounts, Enots)

norm = getNormZ(zPos, zNeg)

ma = np.argsort(norm)


print np.array(kmers)[ma][-10:]
print np.array(norm)[ma][-10:]
