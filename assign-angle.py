import numpy as np
import argparse
import matplotlib.pyplot as plt

# parser = argparse.ArgumentParser(description="")
# parser.add_argument("-bias", type=str, help="bias value to analyse")
# parser.add_argument("-clean", action='store_true', help="whether to clean data or not")
# parser.add_argument("-row", choices=["x","z"], help="which direction to average out")
# parser.add_argument("-remove_NN", action='store_true', help="whether to remove sites with less than 4 NN's")
# args = parser.parse_args()

data = np.genfromtxt('thetasoft.-0.7250.txt', delimiter=' ')
data = data.reshape((-1,20,12))

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

data_txz = np.zeros((data.shape[0]-3000, data.shape[1], data.shape[2]))
data_txz[:, ::2, :] = data[3000:, 0:10, :]
data_txz[:, 1::2, :] = data[3000:, 10:20, :]
ligdata = data_txz
T = len(ligdata)

# get theta values at each lattice site as a function of time
theta_lat = []
for i in range(T):
    theta_lat.append( ligdata[i, :, :].flatten() )
theta_lat = np.array(theta_lat)

# make coarse grained lattice
cd_data = np.genfromtxt('paired-sites.txt')
n_cd = len(cd_data)
# print n_cd
cg_lat = []

for i in range(n_cd):
    cg_lat.append((cd_data[i][2] + 4.12*0.5, cd_data[i][3], cd_data[i][4] + 6.75*0.5))

cg_lat = np.array(cg_lat)
cg_xz = np.transpose( np.array((cg_lat[:, 0], cg_lat[:, 2])) )
# print cg_xz
nn_dist = np.sqrt( (4.12*0.5) **2 + (6.75*0.5)**2 )
# print cd_data[:,2].shape

full_indices = []
indices = []
# now loop over coarse grained lattice and get nearest neighbours for each point
for j in range(n_cd):
    x = cg_xz[j][0]
    z = cg_xz[j][1]
    # now find the Cd atoms closest to this coarse grained site
    dist = np.sqrt( (x - cd_data[:,2])**2 + (z - cd_data[:,4])**2 )
    full_indices.append( np.where( dist <= nn_dist+0.05 ) )
#     print j, indices

# remove_NN = args.remove_NN
remove_NN = None
if remove_NN == None:
    for j in range(n_cd):
        indices.append( full_indices[j] )
else:
    for j in range(n_cd):
        if ( np.array(full_indices[j]).shape[1] == 4 ):
            indices.append( full_indices[j] )        

# now assign theta values to each coarse grained site as a function of time
cg_theta = []
for i in range(T):
    theta_t = []
    for j in range(len(indices)):
        theta_site = theta_lat[i][indices[j]]
        theta_t.append(np.mean(theta_site))
    cg_theta.append(np.array(theta_t))

cg_theta = np.array(cg_theta)
if remove_NN == None:
    cg_theta = cg_theta.reshape((-1,20,12))

bins = np.linspace(-0.9, -0.5, 200)
cg_mean = np.mean(cg_theta, axis=1)
cg_hist, bins = np.histogram(cg_mean, bins=bins)
# cg_hist, bins = np.histogram(cg_theta[1000], bins=bins)
bins = 0.5 * (bins[1:] + bins[:-1])
plt.plot(bins, cg_hist, 'ro')
plt.plot(bins, cg_hist)
plt.show()




