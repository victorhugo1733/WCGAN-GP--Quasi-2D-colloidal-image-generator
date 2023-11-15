from numpy.core.fromnumeric import shape
import pandas as pd
import numpy as np
import os
import time
from multiprocessing import Pool
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt


def Swap(arr, start_index, last_index): 
    arr[:, [start_index, last_index]] = arr[:, [last_index, start_index]]
    
    
num_frame=1

for i in range(0, 1000):
    globals()["list_data_" + str(i)] = []
    globals()["data_" + str(i)]  = pd.read_csv(r'D:\Images_Polimero0 NaCl0\Victor\Victor\Hugo\Hugo\800x800\Generadas\frames\0.15\frame{i}.csv'.format(i=i))
    globals()["data_" + str(i)]=np.array(globals()["data_" + str(i)])
    Swap(globals()["data_" + str(i)], 1, 2) 
    globals()["data_" + str(i)]=np.delete(globals()["data_" + str(i)], 0, axis=1)


def rdf(particles, dr, rho=None, rcutoff=0.2, eps=1e-15, progress=False):
    start = time.time()

    mins = np.array([0, 0]) #obtenemos el minimo en x y y (es una matriz de dos valores)
    maxs = np.array([600, 600]) #obtenemos el maximo en x y y (es una matriz de dos valores)
    # translate particles such that the particle with min coords is at origin
    particles = particles - mins

    # dimensions of box
    dims = maxs - mins
    
    r_max = (np.min(dims) / 2)*rcutoff
    radii = np.arange(dr, r_max, dr)

    N, d = particles.shape
    if not rho:
        rho = N / np.prod(dims) # number density
    
    # create a KDTree for fast nearest-neighbor lookup of particles
    tree = cKDTree(particles)

    g_r = np.zeros(shape=(len(radii)))
    for r_idx, r in enumerate(radii):
        # find all particles that are at least r + dr away from the edges of the box
        valid_idxs = np.bitwise_and.reduce([(particles[:, i]-(r+dr) >= mins[i]) & (particles[:, i]+(r+dr) <= maxs[i]) for i in range(d)])
        valid_particles = particles[valid_idxs]
            
        # compute n_i(r) for valid particles.
        for particle in valid_particles:
            n = tree.query_ball_point(particle, r+dr-eps, return_length=True) - tree.query_ball_point(particle, r, return_length=True)
            g_r[r_idx] += n
            
        # normalize
        n_valid = len(valid_particles)
        shell_vol = (4/3)*np.pi*((r+dr)**3 - r**3) if d == 3 else np.pi*((r+dr)**2 - r**2)
        g_r[r_idx] /= n_valid*shell_vol*rho 

        if progress:
            print('Computing RDF     Radius {}/{}    Time elapsed: {:.3f} s'.format(r_idx+1, len(radii), time.time()-start), end='\r', flush=True)

    return g_r, radii

A=0
g_r_015g_values = []
for iss in range(1000):
    B=A+1
    for i in range(A, B):
        print('frame=',i)
        globals()["g_r_" + str(i)], globals()["radii_" + str(i)] = rdf(globals()["data_" + str(i)], dr=0.5)
        g_r_015g = np.zeros(shape=(len(g_r_0)))    
        for j in range(len(g_r_0)):
            for i in range(A, B):
                g_r_015g [j] +=globals()["g_r_" + str(i)][j]
                g_r_015g [j] /=num_frame
    A=A+1
    g_r_015g_values.append(g_r_015g)
    
  
#------------------------------------------------------
#------------------------------------------------------
#------------------------------------------------------
#------------------------------------------------------
#------------------------------------------------------
g_r_015g = np.zeros(shape=(len(g_r_0)))
for j in range(len(g_r_0)):
    for i in range(num_frame):
        g_r_015g [j] +=globals()["g_r_" + str(i)][j]
    g_r_015g [j] /=num_frame    
#-------------------Grafica

pix_um=2.44042969
r=(radii_0/7)


fig, ax = plt.subplots(figsize=(10,7))
plt.xlim([0, 8])
plt.ylim([0, 2.8])

for iss in range(990):
    plt.plot(r,g_r_005g_values[iss], lw=1, color=(0.7, 0.9, 0.7))#,label=r"$\rho$=1"
for iss in range(990):
    plt.plot(r,g_r_005_values[iss], lw=1,  color=(0.7, 0.7, 0.7))#,label=r"$\rho$=1"    
    

plt.plot(r, g_r_005, lw=3 ,label=r"Experiment")#
plt.plot(0,0, lw=3 ,  color=(0.7, 0.7, 0.7) ,label=r"$\delta_E$")#
plt.plot(r,g_r_005g, lw=3 ,label=r"WCGAN-GP")#
plt.plot(0,0, lw=3 , color=(0.7, 0.9, 0.7),label=r"$\delta_W$")#
plt.title(r"$\rho=0.05$", fontdict={'fontsize': 24, 'fontweight': 'bold'})
plt.ylabel('g (r)', fontdict = {'fontsize':20, 'fontweight':'bold'},)
plt.xlabel('r/$\sigma$', fontdict = {'fontsize':20, 'fontweight':'bold'})
ax.tick_params(axis='both', which='major', labelsize=15)
plot_properties = dict(linewidth=2 )
plt. rc('font', size = 15)
#numeros
plt.tick_params(axis='both', which='major', labelsize=15, width=3,)
plt.rcParams["font.weight"] = "bold"
#plt.title("$\rho$=0.15", fontdict = {'fontsize':36})

# modifico el spines
ax.spines['top'].set_linewidth(3)
ax.spines['right'].set_linewidth(3)
ax.spines['left'].set_linewidth(3)
ax.spines['bottom'].set_linewidth(3)


plt.legend(loc='lower right')
    