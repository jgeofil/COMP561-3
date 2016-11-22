import csv
import numpy as np

o = []
ii = []
m = []

with open('O.csv', 'rb') as csvfile:
     lines = csv.reader(csvfile, delimiter='\t')
     for i, row in enumerate(lines):
        o.append(int(row[0]))


with open('I.csv', 'rb') as csvfile:
     lines = csv.reader(csvfile, delimiter='\t')
     for i, row in enumerate(lines):
        ii.append(int(row[0]))



with open('M.csv', 'rb') as csvfile:
     lines = csv.reader(csvfile, delimiter='\t')
     for i, row in enumerate(lines):
        m.append(int(row[0]))


import matplotlib.pyplot as plt

bins = np.linspace(0, 200, 100)

plt.hist(o, bins)  # plt.hist passes it's arguments to np.histogram
plt.title("Length frequency distribution of hydrophobic regions")
plt.xlabel('Length')
plt.ylabel('Frequency')
plt.show()


plt.hist(ii, bins)  # plt.hist passes it's arguments to np.histogram
plt.title("Length frequency distribution of hydrophilic regions")
plt.xlabel('Length')
plt.ylabel('Frequency')
plt.show()


plt.hist(m, bins)  # plt.hist passes it's arguments to np.histogram
plt.title("Length frequency distribution of mixed regions")
plt.xlabel('Length')
plt.ylabel('Frequency')
plt.show()
