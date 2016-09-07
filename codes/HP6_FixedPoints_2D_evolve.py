# this script is to analyze the fixed points of HNNP and HN6

import HNRG
import numpy as np
import matplotlib.pyplot as plt



def main():
  nt = 5
  t_min = 1.1
  t_max = 100.
  ms = np.linspace(t_min, t_max, num = nt)
  
  # HNNP: fixed point analysis 
  hp = HNRG.hp()              # class object of HNNP
  m = 3.1
  n = 1000000
  l = n - 10000
  nb = 30
  nls = 4
  ks = hp.rg_n_k(m, n, l)
  plt.figure(1, figsize = (8, 8))
  for i in range(nls):
    hist = HNRG.hist_fx(ks[i*l/nls:(i+1)*l/nls], nbins = nb, density = True)
    plt.subplot(int(np.sqrt(nls)), int(nls/np.sqrt(nls)), 1)
    plt.plot(hist[:,0], hist[:, 1], '-', linewidth = 2.0, label = str(i*l/nls+n-l) + '~' + str((i+1)*l/nls+n-l))
    if i>0:
      plt.subplot(int(np.sqrt(nls)), int(nls/np.sqrt(nls)), i+1)
      plt.plot(hist[:,0], hist[:, 1], '-o')
  
  for i in range(nls):
    plt.subplot(int(np.sqrt(nls)), int(nls/np.sqrt(nls)), i+1)
    if i==0:
      plt.xlim([-0.02, 0.4])
    #plt.ylim([np.amin(hist[:,1]*0.96), np.amax(hist[:,1])*1.04])
    plt.grid('on')
    #if i == 0:
    #  plt.legend()
    if i%int(nls/np.sqrt(nls))==0:
      plt.ylabel('Density')
    if (i+1) >(int(np.sqrt(nls))-1)*int(nls/np.sqrt(nls)):
      plt.xlabel('values of $\kappa$', fontsize = 16)
    if i>0:
      plt.title('HNNP, RG steps: ' + str(i*l/nls+n-l) + '~' + str((i+1)*l/nls+n-l))
    else:
      plt.title('HNNP, RG steps: ' + str(n-l) + '~' + str(n))

  plt.savefig('HNNP_m'+str(m)+'Nmax' + str(n) + 'Nmin' + str(n-l) + '.pdf')

  # HN6: fixed point analysis
  h6 = HNRG.h6()              # class object of HN6
  
  return 
  


  
if __name__ == "__main__":
  main()

