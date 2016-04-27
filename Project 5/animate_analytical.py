import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation



N = 1000
x = np.linspace(0, 1, N)

u_s = 1-x

  
def u_ana(x,t):
  M = 1000
  u = np.zeros(N)
  for k in range(1,M):
    u += -2/(k*np.pi)*np.sin(k*x*np.pi)*np.exp(-(np.pi*k)*t)
  return u


fig, ax = plt.subplots()

line, = ax.plot(x, u_s, linewidth=2)
plt.xlabel('x')
plt.ylabel('Concentration')

plt.ylim(-1, 1)
plt.xlim(0, 1)

def animate(t):
  line.set_ydata(u_ana(x,t)+u_s)
  return line,


ani = animation.FuncAnimation(fig, animate, np.linspace(0, 2, 300), interval=100)

plt.show()