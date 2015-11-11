import matplotlib.pyplot as plt
from numpy import *

data = loadtxt('therm_data20.txt',skiprows=0)
T = array(data[1:,0])
C_v20 = array(data[1:,3])

data = loadtxt('therm_data40.txt',skiprows=0)
C_v40 = array(data[1:,3])

data = loadtxt('therm_data60.txt',skiprows=0)
C_v60 = array(data[1:,3])

data = loadtxt('therm_data80.txt',skiprows=0)
C_v80 = array(data[1:,3])

N = int(array(data[0,1]))

exact = 2.269

plt.plot(T,C_v20, label='$C_v, L=20$')
plt.plot(T,C_v40, label='$C_v, L=40$')
plt.plot(T,C_v60, label='$C_v, L=60$')
plt.plot(T,C_v80, label='$C_v, L=80$')
plt.plot((exact, exact), (0, 2.5),'m--' , label='$J/T_C k_B$' )
plt.title('$C_v (T)$ per particle for different lattice sizes. No. cycles: {}'.format(N), fontsize=18)
plt.xlabel('$T`$ [$J/ T$ $k_B$]')
plt.axis([2, 2.7, 0, 2.5]) 
plt.legend()
plt.show()
