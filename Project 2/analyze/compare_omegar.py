import subprocess
import os
import numpy as np
import matplotlib.pyplot as plt


text_file = 'test.txt'
executable = os.path.join('..','./project2')
target_folder = os.path.join('..','results')

def get_E():
  with open(text_file) as f:
    lines = f.readlines()
  
  E = [0]*3
  for i in range(3):
    E[i] = float(lines[i+5].strip())
  return E

def get_eigv():
  with open(text_file) as f:
    lines = f.readlines()
  N = int(lines[0].split()[0])
  eigv = [[0]*N,[0]*N,[0]*N];
  for i, j in enumerate(range(9,12)):
    eigv[i] = list(map(float, lines[j].split()))
  return eigv

def eigv_to_psi(eigv):
  return [[x*x for x in eigvi] for eigvi in eigv]



def get_psi(N, rho_inf, omega):
  FNULL = open(os.devnull, 'w')
  command = [executable, str(N), str(rho_inf), '1', str(omega)]
  subprocess.call(command, stdout=FNULL)
  eigv = get_eigv()
  return eigv_to_psi(eigv)




def plot_psi_omegar(omega):
  fig = [plt.figure(), plt.figure(), plt.figure()]
  ax = [0]*3
  for i in range(3):
    ax[i] = fig[i].add_subplot(111)
  
  N = 100
  rho_inf = 6
  rho = np.linspace(0, rho_inf, N)
  
  for omegar in omega:
    psi = get_psi(N, rho_inf, omegar)
    E = get_E()
    for i in range(3):
      ax[i].plot(rho, psi[i], 
                 label='$\omega_r = {}$, $E = {}$'.format(omegar, E[i]))
      
  for i in range(3):
    ax[i].set_title('$\psi_{}$'.format(i))
    ax[i].legend(loc='upper right')
    ax[i].set_xlabel(r'$\rho$')
    ax[i].set_ylabel(r'$|\psi|^2$')
    
    ax[i].set_xlim((0, rho_inf))
    savename = os.path.join(target_folder, 
                            'psi_compare_omegar{}.png'.format(i))
    fig[i].savefig(savename, dpi = 400,
                   bbox_inches='tight')
  
  
omega_r = [0.01, 0.5, 1, 5]

plot_psi_omegar(omega_r)



