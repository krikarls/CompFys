import subprocess
import matplotlib.pyplot as plt
import numpy as np
import os

saveloc = os.path.join('report', 'pics')

def run_command(cmd):
  p = subprocess.Popen(cmd)
  p.wait()

def get_lines():
  with open('test.txt') as f:
    return f.readlines()

def plot_stuff(u_ANA, u_eulerF, u_eulerB, u_CN, u_MC, u_MCG):
  fig, ax = plt.subplots()

  N = len(u_ANA)
  x = np.linspace(0, 1, N)
  plt.plot(x, u_ANA, linewidth=2, label='Analytical solution')
  plt.plot(x, u_eulerF, linewidth=2, label='Euler forward')
  plt.plot(x, u_eulerB, linewidth=2, label='Euler backward')
  plt.plot(x, u_CN, linewidth=2, label='Crank Nicolson')
  plt.plot(x, u_MC, linewidth=2, label='Monte Carlo')
  plt.plot(x, u_MCG, linewidth=2, label='Monte Carlo variable')

  plt.title('Diffusion of particles',fontsize=16)
  plt.xlabel('x')
  plt.ylabel('Concentration')
  ax.axvspan(-0.2, 0, ymin=0, ymax=1, alpha=0.5, color='red')
  ax.axvspan(1, 1.2, ymin=0, ymax=1, alpha=0.5, color='red')

  ax.text(-0.1, 0.6, r'axon (presynaptic)', rotation='vertical', fontsize=14)
  ax.text(1.1, 0.65, r'dendrite (postsynaptic)', rotation='vertical', fontsize=14)

  plt.axis([-0.2, 1.2, 0, 1])

  plt.legend(loc='best')
  savename = os.path.join(saveloc, 'dx={}_t={}.png'.format(dx, T))
  plt.savefig(savename, dpi=400, bbox_inches='tight')
  plt.show()
  
  
alpha = 1.0/2 # Stability criteria
dx = 1/100
dt = dx*dx*alpha
T = 0.2

cmd = ['./project5', str(dt), str(dx), str(T)]

run_command(cmd)
lines = get_lines()



u_ANA = list(map(float, lines[7][1:-2].split()))
u_eulerF = list(map(float, lines[9][1:-2].split()))
u_eulerB = list(map(float, lines[12][1:-2].split()))
u_CN = list(map(float, lines[15][1:-2].split()))
u_MC = list(map(float, lines[18][1:-2].split()))
u_MCG = list(map(float, lines[23][1:-2].split()))

MC_variance = list(map(float, lines[7][1:-2].split()))
MCG_variance = list(map(float, lines[7][1:-2].split()))

E = [float(lines[10][0]), float(lines[13][0]), float(lines[16][0]),
     float(lines[19][0]), float(lines[24][0])]

plot_stuff(u_ANA, u_eulerF, u_eulerB, u_CN, u_MC, u_MCG)
