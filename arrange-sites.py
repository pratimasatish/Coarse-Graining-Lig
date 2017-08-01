import numpy as np

data = np.genfromtxt('surf-Cd.txt')
ligdata = np.genfromtxt('surf-lig.txt')

flip = 0
count = 0
size = len(data)
rows = size/12
for i in range(rows/2):
    start = i*12
#     print start
    for j in range(12):
        print "{:d} {:d} {:2.5f} {:2.5f} {:2.5f} {:d}".format( int(data[start+j][0]),int(data[start+j][1]),data[start+j][2],data[start+j][3],data[start+j][4], int(ligdata[start+j][0]) )
#         print "{:d} {:d} {:2.5f} {:2.5f} {:2.5f}".format( int(ligdata[start+j][0]),int(ligdata[start+j][1]),ligdata[start+j][2],ligdata[start+j][3],ligdata[start+j][4] )
    for j in range(12):
        print "{:d} {:d} {:2.5f} {:2.5f} {:2.5f} {:d}".format( int(data[start+j+size/2][0]), int(data[start+j+size/2][1]),data[start+j+size/2][2],data[start+j+size/2][3],data[start+j+size/2][4], int(ligdata[start+j+size/2][0]) )
#         print "{:d} {:d} {:2.5f} {:2.5f} {:2.5f}".format( int(ligdata[start+j+size/2][0]), int(ligdata[start+j+size/2][1]),ligdata[start+j+size/2][2],ligdata[start+j+size/2][3],ligdata[start+j+size/2][4] )





