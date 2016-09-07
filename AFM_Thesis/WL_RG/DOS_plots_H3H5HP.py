import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rc('text', usetex = True)

def main():
  titlename = 'HNs Anti-Ferr Ising'
  filename = 'HN35P6DOS_aps.png'
  a3 = np.loadtxt("H3AFDOS9_r0")
  a5 = np.loadtxt("H5AFDOS9_r0")
  ap = np.loadtxt("HPAFDOS9_r0")
  a6 = np.loadtxt("H6AFDOS9_r0")
  k=9
  #a8 = np.loadtxt("HPAFDOS8_r0")
  #a9 = np.loadtxt("HPAFDOS9_r0")
  styles = ['--r', '-.b', ':k', '-g']

  plt.figure(1, figsize=(7,5))
  p1=plt.plot(a3[:,0]/2**k,a3[:,1]/2**k,styles[0], markersize=8, linewidth = 3.3)
  p2=plt.plot(a5[:,0]/2**k,a5[:,1]/2**k,styles[1], markersize=8, linewidth = 3.3)
  p3=plt.plot(ap[:,0]/2**k,ap[:,1]/2**k,styles[2], markersize=8, linewidth = 3.3)
  p4=plt.plot(a6[:,0]/2**k,a6[:,1]/2**k,styles[3], markersize=8, linewidth = 3.3)
  plt.ylim((0, np.amax([a5[:,1]/2**k])*1.03))
  #p6=plt.plot(a8[:,0]/2**8,a8[:,1]/2**8,styles[5], markersize=8, linewidth = 3.3)
  #p7=plt.plot(a9[:,0]/2**9,a9[:,1]/2**9,styles[6], markersize=8, linewidth = 3.3)
  plt.xlabel('Energy Density $E/N$', fontsize = 18)
  plt.ylabel('$\log(g_E) / N$', fontsize = 18)
  #plt.title(titlename, fontsize = 16)
  plt.grid('on')
  plt.legend([p1[0],p2[0],p3[0], p4[0]],('HN3', 'HN5', 'HNNP', 'HN6'),'upper right')
  plt.xticks(fontsize = 14)
  plt.yticks(fontsize = 14)
  plt.tight_layout()
  #plt.text(-1.4, 0.64, 'N=256', fontsize=24)
  #plt.savefig(filename, dpi=500)
  plt.savefig('HN35P6DOS_h3.pdf')
  

  plt.show()

  
  return 0;
  
  


if __name__ == '__main__':
  main()

