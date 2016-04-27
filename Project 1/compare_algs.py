#!/usr/bin/env python3
import numpy as np
import subprocess
import matplotlib.pylab as plt
import os

filename = 'test.txt'

N = [10, 25, 50, 100, 250, 500, 1000, 2000, 3000, 4000, 5000, 6000]

Alg = ['M','T','S','L']

def analyze_file(filename, N):
  with open(filename) as f:
    for ii in range(2*N+4):
      f.readline()
    time = float(f.readline().split()[0])
    err  = float(f.readline().split()[0])
  return time, err

time = np.zeros((len(Alg),len(N)))
err = np.zeros((len(Alg),len(N)))

for i, N_ in enumerate(N):
	for j, A in enumerate(Alg):
		command = ['./project1', str(N_), filename, A]
		subprocess.call(command)
		time[j,i], err[j,i] = analyze_file(filename, N[i])


log_n = np.log10(N)


plt.plot(log_n, err[0,:], label="Gaussian elimination")
plt.plot(log_n, err[1,:], label="Thomas algorithm")
plt.plot(log_n, err[3,:], label="LU decomposition")
plt.plot(log_n, err[2,:], label="Sparse matrix decomposition")
plt.xlabel('Log10 number of messh points')
plt.ylabel('Log10 relative error')
plt.title('Max relative error for different algorithms')
plt.legend()
plt.savefig(os.path.join('Results', 'n_err.png'), dpi=600, bbox_inches='tight')
plt.show()

log_time = np.log10(time)

plt.plot(log_n, log_time[0,:], label="Gaussian elimination")
plt.plot(log_n, log_time[1,:], label="Thomas algorithm")
plt.plot(log_n, log_time[3,:], label="LU decomposition")
plt.plot(log_n, log_time[2,:], label="Sparse matrix decomposition")
plt.xlabel('$\log_{10}$ number of mesh points')
plt.ylabel('$\log_{10}$ Effective seconds [Clock cycles/MHz]')
plt.title('Comparing time')
plt.legend(loc='upper left')
plt.savefig(os.path.join('Results', 'n_time.png'), dpi=600, bbox_inches='tight')
plt.show()


np.save(os.path.join('Results', 'time.npy'), time)

np.save(os.path.join('Results', 'err.npy'), err)

