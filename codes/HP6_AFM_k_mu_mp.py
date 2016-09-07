# Python code: HN3 and HN5 Phase diagram
import numpy as np
import matplotlib.pyplot as plt
import mpmath as mp
import HNRG as hnrg

def main():
  nt = 80
  ms = np.linspace(0.01, 1, num = nt)
  n = 10000
  l = 100
  ys = np.array([0, 0.5, 1, 1.5, 2])
  y = 0
  plt.figure(1, figsize = (6,16))
  hp = hnrg.hp6mp(100)
  #ks = np.zeros((num_init,))
  for yi in range(ys.size):
    y = ys[yi]
    plt.subplot(ys.size, 1, yi+1)
    for i in range(ms.size):
      m =1./ ms[i]
      ks = np.unique(np.array(hp.rg_n_k(mp.mpf(m), mp.mpf(y), n, l), dtype=float))
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
  plt.savefig('HP6_AFM_k_mu_N'+str(n)+'l'+str(l)+'ys_mp_p.pdf')
  plt.show()
  
if __name__ == '__main__':
  main()
  
