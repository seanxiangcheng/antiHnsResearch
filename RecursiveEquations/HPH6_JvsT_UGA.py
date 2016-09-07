# Python code: HN3 and HN5 Phase diagram
import numpy as np
import matplotlib.pyplot as plt
import HNRG as hnrg

def main():
  nt = 120
  ts = np.linspace(0.02, 3.5, num = nt)
  ms = np.exp(2.0/ts)
  n = 1000
  l = 200
  ys = np.array([0, 0])
  y = 0
  plt.figure(1, figsize = (7,7))
  hp = hnrg.hp6()
  #ks = np.zeros((num_init,))
  for yi in range(1):
    y = ys[yi]
    plt.subplot(1, 1, yi+1)
    for i in range(ms.size):
      m =ms[i]
      t = ts[i]
      js = -np.log(np.unique(hp.rg_n_k(m, y, n, l)))*(1.0)*t/4.0
      plt.plot(t*np.ones((js.size,)), js, 'ks', markersize = 0.8)
      if js.size<3:
        plt.plot(t*np.ones((js.size,)), js, 'bs', markersize = 1)
      else:
        plt.plot(t*np.ones((js.size,)), js, 'ks', markersize = 0.8)
      #plt.plot(ms[i]*np.ones((ks.size,)), ks, 'ks', markersize = 1)

  #plt.semilogy(np.linspace(0.3333, 1, num = nt), np.ones((nt,)), 'k-', linewidth = 2)
    plt.grid('on')
    #plt.xlim([0., 1.01])
    #plt.ylim([0, 1.01])
    #plt.legend(loc = 2)
    if yi==0:
      plt.ylabel('$J$', fontsize = 16)
    plt.xlabel('$T$', fontsize = 16)
    #plt.title('HP6,y='+str(y)+',steps:'+str(n-l)+'~'+str(n), fontsize = 16)
  #plt.savefig('HP6_AFM_J_T_ys_np.pdf')
  #plt.savefig('HP6_AFM_J_T_ys_np.png', dpi = 100)
  """
  jfps = np.zeros((nt,))
  t_trans = [1.82, 1.5]
  for yi in range(1):
    y = ys[yi]
    plt.subplot(int(np.ceil(ys.size/2.)), 1, yi+1)
    for i in range(ts.size):
      t = ts[i]
      m = ms[i]
      jfps[i] = -np.log(hp.kappa_fp(m, y))*(1.0)*t/4.0
      stable_index=0
      for stable_index in range(jfps.size):
        if ts[stable_index]> t_trans[yi]:
          break
    plt.plot(ts[stable_index:], jfps[stable_index:], '-k', linewidth = 2.5, label = "RG stable")
    plt.plot(ts[:stable_index+2], jfps[:stable_index+2], '--r', linewidth = 2.5, label = "RG unstable")
    plt.plot([],[],'ks', markersize = 3, label = 'RG Numerical')
    plt.legend(loc='upper right')
  """
  plt.plot([],[],'ks', markersize = 3, label = 'RG Numerical')
  plt.legend(loc='upper right')
  #plt.savefig("HPH6_JvsT_aps1.png", dpi=500)
  #plt.savefig("HPH6_JvsT_aps1.pdf")

  plt.show()
  
if __name__ == '__main__':
  main()
  
