# this script is to analysis the stability from the eigenvalues
import HNRG 
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
def main():
  hp = HNRG.hp6()
  N = 100

  styles = ['k', 'r:', 'y','g:', 'c', 'm:','b', 'k:', 'r']
  eigs_real = np.zeros((N, 2))
  eigs_abs = np.zeros((N, 2))
  jac_det = np.zeros((N,))
  jac_tr = np.zeros((N,))
  ys = np.array([-2., -1.5, -1., -0.5, 0.,0.5, 1.,1.5, 2.])
  ms = np.linspace(0.02, 1, num = N)
  plt.figure(1, figsize = (10,12))

  for j in range(ys.size):
    y = ys[j]
    for i in range(N):
      m = 1./ms[i]
      eigs_real[i,:] = hp.jac_eig_real(m, y)
      eigs_abs[i, :] = hp.jac_eig_abs(m, y)
    plt.subplot(2,2,1)
    plt.plot(ms, eigs_real[:, 0], styles[j], linewidth=2.5, label = 'y='+str(y))
    plt.plot(ms, eigs_real[:, 1], styles[j], linewidth=2.5)
    plt.subplot(2,2,2)
    plt.plot(ms, eigs_abs[:, 0], styles[j], linewidth=2.5, label = 'y='+str(y))
    plt.plot(ms, eigs_abs[:,1],  styles[j], linewidth=2.5)
  plt.subplot(2,2,1)
  plt.grid('on')
  plt.xlabel('1/$\mu$')
  plt.ylabel('Real(eigenvalues)')
  plt.title('Real part of eigenvalues')
  #plt.legend()
  plt.subplot(2,2,2)
  plt.grid('on')
  plt.xlabel('1/$\mu$')
  plt.ylabel('Abs(eigenvalues)')
  plt.title('Magnitude of eigenvalues')
  #plt.legend()
  #plt.savefig("HP6_Eigen_plot.pdf")
  
  for j in range(ys.size):
    y = ys[j]
    for i in range(N):
      m = 1./ms[i]
      jac_mat = hp.jacobian_matrix2(m, y)
      jac_det[i] = np.linalg.det(jac_mat)
      jac_tr[i] = jac_mat[0,0] + jac_mat[1,1]
    plt.subplot(2,2,3)
    plt.plot(ms, jac_det, styles[j], linewidth=2.5, label = 'y='+str(y))
    plt.subplot(2,2,4)
    plt.plot(jac_tr, jac_det, styles[j], linewidth=2.5, label = 'y='+str(y))
  plt.subplot(2,2,3)
  plt.legend()
  plt.xlabel('1/$\mu$')
  plt.ylabel('Determinant of Jacobian')
  plt.title('Determinant of Jacobian')
  plt.grid('on')
  plt.subplot(2,2,4)
  plt.grid('on')
  plt.xlabel('Trace')
  plt.ylabel('determinant of Jacobian')
  plt.title('Det vs Tr')
  plt.legend()
  #plt.savefig("HP6_Eigen_plot2.pdf")
  plt.show()
  return ;


if __name__ =="__main__":
  main()
