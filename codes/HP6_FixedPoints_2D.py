# this script is to analyze the fixed points of HNNP and HN6
"""
1. check consistence over RG
2. 2D plot: hist(k)
3. 2D plot: k vs mu
"""
import HNRG
import numpy as np
import matplotlib.pyplot as plt



def main():
  nt = 4
  t_min = 0.18
  t_max = 0.3
  ys = np.array([0, 0.5, 1, 2])
  ms = np.linspace(t_min, t_max, num = nt)
  nb = 120
  n = 10000
  l = n - 100
  k_mm = np.array([0, 100])
  # styles parameters#
  msize = 4
  lwidth = 2.0
  styles = ['-o', '-->']
  plt.figure(1, figsize=(16,5))
  hp = HNRG.hp6()              # class object of HNNP
  for yi in range(ys.size):
    plt.subplot(1, ys.size, yi+1)
    y = ys[yi]
    for i in range(nt):
      m = 1./ms[i]
      print m
      ks = hp.rg_n_k(m, y, n, l)
      hist = HNRG.hist_fx(ks, nbins = nb)
      if m<3.0:
        s = styles[1]
      else: s = styles[0]
      plt.plot(hist[:,0], hist[:,1], s, markersize = msize, linewidth = lwidth,
        label = '$1/\mu$='+str(round(1.0/m, 2)))
  #plt.xlim([-0.2, 200])
  #plt.ylim([0.0, 0.72])
    plt.grid('on')
    if yi==0:
      plt.legend(loc='upper left')
    plt.title('y='+str(y))
    if yi==0:
      plt.ylabel('Density')
    plt.xlabel('tanh($\kappa$)')
  plt.savefig('HP6_tanhFP_freq_ys_p.pdf')
  plt.show()
  return 
  
  
if __name__ == "__main__":
  main()

