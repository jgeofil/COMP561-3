f =  open('GATA2_chr1.fa')
lines = f.read().split('>')[1:]

seqs = [''.join(l.splitlines()[1:]) for l in lines]

kmers = []

def rec (s):
    if len(s) >= 6:
        kmers.append(s)
    else:
        for c in ['A', 'T', 'G', 'C', 'X', 'Y', 'Z']:
            rec(s + c)

rec('')

print len(kmers)


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


t = Tree()
k = 6

for s in seqs:
    for i in range(len(s)-k):
        t.addCount(s[i:i+k])


counts = []

for k in kmers:
    counts.append(t.getCount(k))

print counts
