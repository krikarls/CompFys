from numpy import *
import matplotlib.pyplot as plt

r1 = linspace(0,5,100)
r2 = linspace(-5,0,100)

wave_function1 = exp(-4*(r1))
wave_function2 = exp(4*(r2))

plt.plot(r1,wave_function1, color = 'b')
plt.plot(r2,wave_function2, color = 'b')
plt.title('Wave function, $\psi$', fontsize=18)
plt.xlabel('$r_1 , r_2$', fontsize=16)
plt.show()