# this script is to analyze the fixed points of HNNP and HN6

import HNRG
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter


def main():
  nt = 50
  t_min = 0.2
  t_max = 0.9
  ms = np.linspace(t_min, t_max, num = nt)
  k_mm = np.array([0, 10])
  nb = int(k_mm[1] / (0.008357 + np.random.uniform()*0.001347))
  n = 50000
  l = 1000
  
  hp = HNRG.hp()              # class object of HNNP
  X = ms
  ks = hp.rg_n_k(4, 20000, 500)
  hist = HNRG.hist_fx(ks, nbins = nb, mm = k_mm)
  Y = hist[:,0]
  Z = np.zeros((nb, nt))
  
  #print hp.rg_mm_k(10.01, 100000, 1000)
  #hist =  HNRG.hist_fx(hp.rg_n_k(2.01, 200000, 1000), nbins = k_mm[1]/0.00832, mm = k_mm)

  for i in range(nt):
    m = 1./ms[i]
    ks = hp.rg_n_k(m, n, l)
    hist = HNRG.hist_fx(ks, nbins = nb, mm = k_mm)
    Z[:,i]=hist[:,1]
  
  X, Y = np.meshgrid(X, Y)  #plt.figure(1, figsize = (8, 8))
  
  fig = plt.figure()
  ax = fig.gca(projection='3d')
  surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm,
        linewidth=0, antialiased=False)
  ax.zaxis.set_major_locator(LinearLocator(10))
  ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

  fig.colorbar(surf, shrink=0.5, aspect=5)

  plt.show()
  return 
  


  
if __name__ == "__main__":
  main()

