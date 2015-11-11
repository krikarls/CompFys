import matplotlib.pyplot as plt
from numpy import *

data = loadtxt('therm_data.txt',skiprows=0)

T = array(data[1:,0])
E = array(data[1:,1])
M = array(data[1:,2])
C_v = array(data[1:,3])
Chi = array(data[1:,4])

dim = int(array(data[0,0]))
N = int(array(data[0,1]))

exact = 2.269

max_Cv = argmax(C_v)

print 'The peak of the specific heat is at T= ', T[max_Cv]

plt.plot(T,E, label='$E(T)$')
plt.plot(T,M, label='$M(T)$')
plt.plot(T,C_v, label='$C_v (T)$')
plt.plot(T,Chi, label='$\chi(T)$')
plt.plot((exact, exact), (-2, 2.5),'m--' , label='$J/T_C k_B$' )
plt.title('Thermodynamic properties per particle. Lattice size: {}, cycles: {}'.format(dim, N))
plt.xlabel('$T`$ [$J/ T$ $k_B$]')
plt.legend()
plt.show()
