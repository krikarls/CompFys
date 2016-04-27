#!/usr/bin/env python3

import subprocess

def get_state_size():
  N = 0
  print('Please state the size of the system (> 1)')
  while N<2:
    N = int(input('N = '))
  return str(N)

def get_filename():
  print('Please state the filename to save the file to, default=test.txt')
  filename = input('Filename = ')
  if filename == '':
    filename = 'test.txt'
  return filename

def get_algtype():
  print('Please choose the algorithm to use, M - matrix decomp,\n' +
      'S - sparse matrix decomp, T - Thomas algorithm,\n' +
      'L - LU decomp', 'default=T')
  alg = input('<type> = ')
  if alg not in ['S','M','T', 'L']:
    alg = 'T'
  return alg
  

N = get_state_size()

filename = get_filename()
alg_type = get_algtype()

command = ['./project1', N, filename, alg_type]
subprocess.call(command)

command = ['./plot.py', filename]
subprocess.call(command)
