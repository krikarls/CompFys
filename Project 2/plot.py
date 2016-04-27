import numpy as np
import matplotlib.pyplot as plt

f = open('test.txt')


lines = f.readlines()

N = int(lines[0].split()[0])

rhos = lines[1].split()
rho_0 = float(rhos[0])
rho_inf = float(rhos[1])

rho = np.linspace(rho_0, rho_inf, N)

two_electrons = (lines[2].split()[0] == 'true')
omega_r = float(lines[3].split()[0])

E = [float(lines[i].split()[0]) for i in range(5, 8)]


plt.xlim(rho_0, rho_inf)
plt.xlabel(r'$\rho$')
plt.ylabel('Probability density $|\psi|^2$')

if two_electrons:
  plt.title('Two electrons with coupling $\omega_r = {}$'.format(omega_r))
else:
  plt.title('Radial part of Schrödinger equation in a harmonic oscillator')

psi = [[],[],[]];
for i in range(3):
  psi[i] = [float(x)**2 for x in lines[9+i].split()]
  plt.plot(rho, psi[i], linewidth=2,
	   label='$\psi_{0}, E_{0} = {1}$'.format(i, E[i]))
  

plt.legend(loc='upper right')
plt.show()


