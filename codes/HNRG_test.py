#module test
import HNRG 
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt

def main():
  
  
  hp = HNRG.hp6()
  """
  # This block is to plot k vs m
  N = 50
  ys = np.array([-2., -1., 0., 1., 2.])
  kaps = np.zeros((N, ys.size))
  ms = np.linspace(0.02, 1, num = N)
  plt.figure(1, figsize=(6,6))
  for j in range(ys.size):
    y = ys[j]
    for i in range(N):
      m = 1./ms[i]
      kaps[i, j] = hp.kappa_fp(m, y)
    plt.plot(ms, kaps[:,j],'-')
  plt.ylim([0, 6.5])
  plt.show()
  """
  N = 50
  y = 0.0
  eigs_real = np.zeros((N, 2))
  eigs_abs = np.zeros((N, 2))
  ms = np.linspace(0.02, 1, num = N)
  for i in range(N):
    m = 1./ms[i]
    eigs_real[i,:] = hp.jac_eig_real(m, y)
    eigs_abs[i, :] = hp.jac_eig_abs(m, y)
  
  plt.figure(1, figsize = (6,12))
  plt.subplot(2,1,1)
  plt.plot(ms, eigs_real[:, 0], linewidth=3)
  plt.subplot(2,1,2)
  plt.plot(ms, eigs_abs[:, 0], linewidth=3)
  plt.show()
  return ;
  
  
if __name__ == "__main__":
  main()
  
