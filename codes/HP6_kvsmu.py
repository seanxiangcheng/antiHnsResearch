# Python code: HN3 and HN5 Phase diagram
import numpy as np
import matplotlib.pyplot as plt
import HNRG as hnrg

def main():
  nt = 40
  ms = np.linspace(0.01, 0.6, num = nt)
  n = 10000000
  l = 200
  ys = np.array([-1, -0.5, 0, 0.5, 1])
  y = 0
  plt.figure(1, figsize = (6,16))
  hp = hnrg.hp6()
  #ks = np.zeros((num_init,))
  for yi in range(ys.size):
    y = ys[yi]
    plt.subplot(ys.size, 1, yi+1)
    for i in range(ms.size):
      m =1/ ms[i]
      ks = np.unique(hp.rg_n_k(m, y, n, l))
      plt.semilogy(ms[i]*np.ones((ks.size,)), ks, 'ks', markersize = 1)
      #plt.plot(ms[i]*np.ones((ks.size,)), ks, 'ks', markersize = 1)

  #plt.semilogy(np.linspace(0.3333, 1, num = nt), np.ones((nt,)), 'k-', linewidth = 2)
    plt.grid('on')
    #plt.xlim([0., 1.01])
    #plt.ylim([0, 1.01])
    #plt.legend(loc = 2)
    plt.ylabel('$\kappa$', fontsize = 16)
    if yi==ys.size-1:
      plt.xlabel('$1/\mu$', fontsize = 16)
    plt.title('HP6,y='+str(y)+',steps:'+str(n-l)+'~'+str(n), fontsize = 16)
  plt.savefig('HP6_AFM_k_mu_N'+str(n)+'l'+str(l)+'ys_big_n3.pdf')
  plt.show()
  
if __name__ == '__main__':
  main()
  
