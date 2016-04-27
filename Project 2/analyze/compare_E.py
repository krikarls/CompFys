import subprocess
import os
import numpy as np
import matplotlib.pyplot as plt


text_file = 'test.txt'
executable = os.path.join('..','./project2')
target_folder = os.path.join('..','results')

def get_E(filename=text_file):
  with open(filename) as f:
    lines = f.readlines()
  
  E = [0]*3
  for i in range(3):
    E[i] = float(lines[i+5].strip())
  return E



def get_E_matrix(Ns, rho_infs):
  FNULL = open(os.devnull, 'w')
  E = np.zeros((len(rho_infs),len(Ns), 3))
  for i, rho_inf in enumerate(rho_infs):
    for j, N in enumerate(Ns):
      command = [executable, str(N), str(rho_inf), str(0)]
      subprocess.call(command, stdout=FNULL)
      tempE = get_E()
      for k in range(3):
        E[i,j,k] = tempE[k]
  return E
      
def plt_rel_error(E, Ns, rho_infs, index):
  numbers = ['first', 'second', 'third']
  E_ex = [3, 7, 11] # Exact solutions
  relE = abs(E[:,:,index]-E_ex[index])/E_ex[index]

  for i, rho_inf in enumerate(rho_infs):
    plt.semilogy(Ns, relE[i], label=r'$\rho_\infty = {}$'.format(rho_inf),
	     linewidth=2)
    
  plt.legend(loc='best')
  plt.ylim([0,0.5])
  plt.xlim([min(Ns), max(Ns)])
  plt.xlabel('Number of mesh points')
  plt.ylabel('Relative error of {} eigenvalue'.format(numbers[index]))
  plot_name = os.path.join(target_folder,'rel_logE{}.png'.format(index))
  plt.savefig(plot_name , dpi=400, bbox_inches='tight')
  plt.show()
  

rho_infs = [2, 3, 4, 5, 6, 8]
Ns = [10, 25, 50, 75, 100, 125, 150, 200]
E = get_E_matrix(Ns, rho_infs)

for i in range(3):
  plt_rel_error(E, Ns, rho_infs, i)
  
  