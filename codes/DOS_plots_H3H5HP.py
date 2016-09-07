import sys
import numpy as np
import matplotlib.pyplot as plt

def main():
  titlename = 'HNs Anti-Ferr Ising'
  filename = 'HN35P6DOS_aps.png'
  a3 = np.loadtxt("H3AFDOS8_r0")
  a5 = np.loadtxt("H5AFDOS8_r0")
  ap = np.loadtxt("HPAFDOS8_r0")
  a6 = np.loadtxt("H6AFDOS8_r0")
  #a8 = np.loadtxt("HPAFDOS8_r0")
  #a9 = np.loadtxt("HPAFDOS9_r0")
  styles = ['--b', '--r', '-k', '-g']

  plt.figure(1, figsize=(7,5))
  p1=plt.plot(a3[:,0]/2**8,a3[:,1]/2**8,styles[0], markersize=8, linewidth = 3.3)
  p2=plt.plot(a5[:,0]/2**8,a5[:,1]/2**8,styles[1], markersize=8, linewidth = 3.3)
  p3=plt.plot(ap[:,0]/2**8,ap[:,1]/2**8,styles[2], markersize=8, linewidth = 3.3)
  p4=plt.plot(a6[:,0]/2**8,a6[:,1]/2**8,styles[3], markersize=8, linewidth = 3.3)
  plt.ylim((0, np.amax([a5[:,1]/2**8])*1.03))
  #p6=plt.plot(a8[:,0]/2**8,a8[:,1]/2**8,styles[5], markersize=8, linewidth = 3.3)
  #p7=plt.plot(a9[:,0]/2**9,a9[:,1]/2**9,styles[6], markersize=8, linewidth = 3.3)
  plt.xlabel('Energy Density $E/N$', fontsize = 18)
  plt.ylabel('$\log(g) / N$', fontsize = 18)
  #plt.title(titlename, fontsize = 16)
  plt.grid('on')
  plt.legend([p1[0],p2[0],p3[0], p4[0]],('HN3', 'HN5', 'HNNP', 'HN6'),'upper right')
  plt.text(-1.4, 0.64, 'N=256', fontsize=24)
  #plt.savefig(filename, dpi=500)
  #plt.savefig('HN35P6DOS_h3.pdf')
  

  plt.show()

  
  return 0;
  
  


if __name__ == '__main__':
  main()

