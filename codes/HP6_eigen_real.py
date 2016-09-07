# this script is to analysis the stability from the eigenvalues
import HNRG 
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
def main():
  hp = HNRG.hp6()
  N = 100

  styles = ['g', 'r--', 'y','g--', 'c', 'm--','b', 'c--', 'r']
  eigs_highT = np.zeros((N, 2))
  eigs_abs = np.zeros((N, 2))
  jac_det = np.zeros((N,))
  jac_tr = np.zeros((N,))
  ys = np.array([-2., -1.5, -1., -0.5, 0.,0.5, 1.,1.5, 2.])
  ms = np.linspace(0.02, 1, num = N)
  plt.figure(1, figsize = (6,6))

  for j in range(ys.size):
    y = ys[j]
    for i in range(N):
      m = 1./ms[i]
      eigs_abs[i, :] = hp.jac_eig_abs(m, y)
      eigs_highT[i,:] = hp.jac_eig_highT(m, y)
    plt.plot(ms, eigs_abs[:, 0], styles[j], linewidth=2.5, label = 'y='+str(y))
    plt.plot(ms, eigs_abs[:,1],  styles[j], linewidth=2.5)
    plt.plot(ms, eigs_highT[:,0], 'k', linewidth = 2.5)
    #plt.plot(ms, eigs_highT[:,1], 'k', linewidth = 2.5)

  print eigs_highT
  plt.grid('on')
  plt.xlabel('1/$\mu$')
  plt.ylabel('Abs(eigenvalues)')
  plt.title('Magnitude of eigenvalues')
  #plt.legend()
  #plt.savefig("HP6_Eigen_plot.pdf")
  plt.show()
  return ;


if __name__ =="__main__":
  main()
