import numpy as np

data = np.genfromtxt('first-frame.txt')

for i in range(3841, 6001, 18):
    print "{:d} {:d} {:2.6f} {:2.6f} {:2.6f}".format(i, int(data[i][0]), data[i][1], data[i][2], data[i][3])

for i in range(8161, 10321, 18):
    print "{:d} {:d} {:2.6f} {:2.6f} {:2.6f}".format(i, int(data[i][0]), data[i][1], data[i][2], data[i][3])

