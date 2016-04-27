import subprocess
import os
import numpy as np
import matplotlib.pyplot as plt


text_file = 'test.txt'
executable = os.path.join('..','./project2')
target_folder = os.path.join('..','results')


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



def get_psi(N, rho_inf):
  FNULL = open(os.devnull, 'w')
  command = [executable, str(N), str(rho_inf), '1']
  subprocess.call(command, stdout=FNULL)
  eigv = get_eigv()
  return eigv_to_psi(eigv)




def plot_psi_N(Ns):

  rho_inf = 5
  fig = [plt.figure(), plt.figure(), plt.figure()]
  ax = [0]*3
  for i in range(len(ax)):
    ax[i] = fig[i].add_subplot(111)
    
  for N in Ns:
    rho = np.linspace(0, rho_inf, N)
    psi = get_psi(N, rho_inf)
    for i in range(3):
      ax[i].plot(rho, psi[i], linewidth=2,
               label='$N = {}$'.format(N))
  
  for i in range(3):
    ax[i].set_title('$\psi_{}$'.format(i))
    ax[i].legend(loc='upper right')
    ax[i].set_ylabel(r'$|\psi|^2$')
    ax[i].set_xlabel(r'$\rho$')
    savename = os.path.join(target_folder, 'N_compare_psi{}.png'.format(i))
    fig[i].savefig(savename, dpi = 400,
                   bbox_inches='tight')

def plot_psi_rho(rhos):
  N = 100
  fig = [plt.figure(), plt.figure(), plt.figure()]
  ax = [0]*3
  for i in range(len(ax)):
    ax[i] = fig[i].add_subplot(111)
    
  for rho_inf in rhos:
    rho = np.linspace(0, rho_inf, N)
    psi = get_psi(N, rho_inf)
    for i in range(3):
      ax[i].plot(rho, psi[i], linewidth=2,
               label=r'$\rho_\infty = {}$'.format(rho_inf))
      
  for i in range(3):
    ax[i].set_title('$\psi_{}$'.format(i))
    ax[i].legend(loc='upper right')
    ax[i].set_ylabel(r'$|\psi|^2$')
    ax[i].set_xlabel(r'$\rho$')
    ax[i].set_xlim((0,max(rhos)))
    savename = os.path.join(target_folder, 
                            'psi_inf_compare_psi{}.png'.format(i))
    fig[i].savefig(savename, dpi = 400,
                   bbox_inches='tight')
  
  
Ns = [10, 25, 50, 75, 100]
plot_psi_N(Ns)

rhos = [2, 3, 4, 5, 6]
plot_psi_rho(rhos)




