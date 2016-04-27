import matplotlib.pyplot as plt
import subprocess
import os
from math import exp, cosh, sinh
import numpy as np

"""
This program is used to recreate the figures from the report. Simply write True instead
of False in the if-loops to create the figures associated with that exercise
"""


filename = 'test.txt'
saveloc = 'pics'  #os.path.join('report','pics')
if (not os.path.exists(saveloc)):
  os.mkdir(saveloc)

def normalise(array, factor):
  for i in range(len(array)):
    array[i] /= factor

def get_stdout(cmd):
  p = subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines=True)
  return p.stdout.read()



def Z(beta):
  return 12 + 4*cosh(8*beta)

def E(beta, Z):
  return -32*sinh(8*beta)/Z

def M(beta, Z):
  return (4*exp(8*beta) + 2)/Z

def Cv(beta, Z):
  return beta**2*(256*cosh(8*beta)/Z - 1024*sinh(8*beta)**2/Z**2)

def chi(beta, Z):
  return beta**2*(-32*(4*exp(8*beta) + 2)*sinh(8*beta)/(4*cosh(8*beta) + 12)**2 + 32*exp(8*beta)/(4*cosh(8*beta) + 12))



#a) + b) 2x2 lattice
if (False):
  N = 50
  Ts = np.linspace(1, 5, N)
  E_measured = [0]*N
  sigmaE = [0]*N
  M_measured = [0]*N
  sigmaM = [0]*N

  for i, T in enumerate(Ts):
    cmd = ['./project4', 'b', str(T)]
    line = get_stdout(cmd)
    E_measured[i], M_measured[i], sigmaE[i], sigmaM[i] = line.split()


  beta = [1/T for T in Ts]
  Z = [Z(beta_) for beta_ in beta]
  cv = [Cv(beta_, Z_)/4 for beta_, Z_ in zip(beta, Z)]
  Chi = [chi(beta_, Z_)/4 for beta_, Z_ in zip(beta, Z)]


  plt.plot(Ts, E_measured,linewidth=2, label='E measured from Ising model')
  E_analytical = [E(beta_, Z_)/4 for beta_, Z_ in zip(beta, Z)]
  plt.plot(Ts, E_analytical,linewidth=2, label='E analytically solved')
  plt.legend(loc='upper left')
  plt.xlabel('Temperature')
  plt.ylabel('Energy per particle')
  pic_filename = os.path.join(saveloc, 'aE.png')
  plt.savefig(pic_filename,dpi=400,bbox_inches='tight')
  plt.show()

  plt.plot(Ts, M_measured,linewidth=2, label='M measured from Ising model')
  M_analytical = [M(beta_, Z_)/2 for beta_, Z_ in zip(beta, Z)]
  plt.plot(Ts, M_analytical,linewidth=2, label='M analytically solved')
  plt.legend(loc='lower left')
  plt.xlabel('Temperature')
  plt.ylabel('Magnetisation per particle')
  pic_filename = os.path.join(saveloc, 'aM.png')
  plt.savefig(pic_filename,dpi=400,bbox_inches='tight')
  plt.show()

  plt.plot(Ts, sigmaE, linewidth=2, label='$C_v$ measured from Ising model')
  plt.plot(Ts, cv, linewidth=2, label='$C_v$ analytically solved')
  plt.xlabel('Temperature')
  plt.ylabel('$C_v$ per particle')
  plt.legend(loc='lower right')
  pic_filename = os.path.join(saveloc, 'acv.png')
  plt.savefig(pic_filename,dpi=400,bbox_inches='tight')
  plt.show()

  plt.plot(Ts, Chi,linewidth=2,label='$\chi$ measured from Ising model')
  plt.plot(Ts, sigmaM,linewidth=2,label='$\chi$ analytically solved')
  plt.xlabel('Temperature')
  plt.ylabel('$\chi$ per particle')
  plt.legend(loc='upper right')
  pic_filename = os.path.join(saveloc, 'achi.png')
  plt.savefig(pic_filename,dpi=400,bbox_inches='tight')
  plt.show()


#c) Correlation time
if (False):
  dim = 20
  for T in [1, 2.4]:
    command = ['./project4', 'c', str(T)]
    subprocess.call(command)
    with open(filename) as f:
      Eo = list(map(float, f.readline().split()))
      Mo = list(map(float, f.readline().split()))
      ao = list(map(float, f.readline().split()))
      Er = list(map(float, f.readline().split()))
      Mr = list(map(float, f.readline().split()))
      ar = list(map(float, f.readline().split()))

    N = len(Eo)

    spins = dim*dim

    normalise(Eo, spins)

    if Mo[-1]< 0:
      normalise(Mo, -spins)
    else:
      normalise(Mo, spins)


    normalise(Er, spins)
    if Mr[-1]<0:
      normalise(Mr, -spins)
    else:
      normalise(Mr, spins)
      
      
    plt.ylim((-2.1,2.1))
    plt.xlim((0,N))

    plt.plot(Eo,label='Energy, ordered start',linewidth=2)
    plt.plot(Er,label='Energy, random start',linewidth=2)

    plt.plot(Mo,'-',label='Magnetism, ordered start',linewidth=2)
    plt.plot(Mr,'-',label='Magnetism, random start',linewidth=2)

    plt.legend(loc='best')
    plt.xlabel('Number of MC cycles')
    plt.ylabel('Property per particle')
    plt.title('Energy and magnetism for ordered/unordered start with lattice size 20x20')
    pic_filename = os.path.join(saveloc, 'cEM{}.png'.format(T))
    plt.savefig(pic_filename,dpi=400,bbox_inches='tight')
    plt.show()

    plt.plot(ao,label='Ordered start',linewidth=2)
    plt.plot(ar,label='Random start',linewidth=2)
    plt.ylim(0, max(ao[-1],ar[-1])*1.1)
    plt.legend(loc='center right')
    plt.xlim((0,N))
    plt.xlabel('Number of MC cycles')
    plt.ylabel('Accumulated accepted flips')
    plt.title('Number of acceptances for ordered/unordered start with lattice size 20x20')
    pic_filename = os.path.join(saveloc, 'cA{}.png'.format(T))
    plt.savefig(pic_filename,dpi=400,bbox_inches='tight')
    plt.show()


#d) Probability for E
if (False):
  for T in [1, 2.4]:
    M = 10000
    cmd = ['./project4', 'd', str(T), str(M)]

    Epair = get_stdout(cmd).split()

    E = [int(pair.split(',')[0]) for pair in Epair]
    p = [int(pair.split(',')[1]) for pair in Epair]

    normalise(p, sum(p))

    plt.xlim(E[0],E[-1])
    plt.plot(E, p,'-o', linewidth=2)
    plt.title('Probability for energy in a 20x20 system with T = {}'.format(T))
    plt.xlabel('Energy')
    plt.ylabel('p(E)')
    pic_filename = os.path.join(saveloc, 'dT={}.png'.format(T))
    plt.savefig(pic_filename, dpi=400, bbox_inches='tight')
    plt.show()

#e) Critical temperature
if (False):  
  for dim in [20, 40, 60, 80, 100, 200]:
    command = ['./project4', 'e', str(dim)]

    subprocess.call(command)

    with open(filename) as f:
      N = int(f.readline().split()[0]);
      f.readline()
      f.readline()
      f.readline()
      f.readline()
      f.readline()
      f.readline()
      f.readline()
      f.readline()
      f.readline() # Input types
      
      T = [0]*N
      E = [0]*N
      M = [0]*N
      Cv = [0]*N
      chi = [0]*N
      ar = [0]*N
      
      for i in range(N):
        line = f.readline();
        T[i], E[i], M[i], Cv[i], chi[i], ar[i] = map(float, line.split())

    f, ax = plt.subplots(2,2,sharex=True)

    ax[0, 0].plot(T, E, label='Energy', linewidth=2)
    ax[0, 0].set_title('$<E>$')

    ax[0, 1].plot(T, M, label='Magnetism', linewidth=2)
    ax[0, 1].set_title(r'$<M>$')

    ax[1, 0].plot(T, Cv, label=r'$C_v$', linewidth=2)
    ax[1, 0].set_title(r'$C_v$')

    ax[1, 1].plot(T, chi, label=r'$\chi$', linewidth=2)
    ax[1, 1].set_title(r'$\chi$')

    for tick in ax[1,0].get_xticklabels():
      tick.set_rotation(45)
      
    for tick in ax[1,1].get_xticklabels():
      tick.set_rotation(45)

    pic_filename = os.path.join(saveloc, 'e{}.png'.format(dim))
    plt.savefig(pic_filename,dpi=400,bbox_inches='tight')
    plt.show()

    