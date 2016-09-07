
import sys
import numpy as np
import matplotlib.pyplot as plt

def main():
  titlename = 'Ground State Anti-Ferr Ising'
  filename = 'LowestEnergy.pdf'
  a = np.loadtxt("Lowest_Energy_SAWL")

  #a8 = np.loadtxt("HPAFDOS8_r0")
  #a9 = np.loadtxt("HPAFDOS9_r0")
  styles = [ '-bo', '-cs', '-kd']
  

  plt.figure(1, figsize=(8,8))
  p1=plt.plot(a[:,0],a[:,1]/2**np.arange(3,11),styles[0], markersize=10, linewidth = 3.3)
  p2=plt.plot(a[:,0],a[:,2]/2**np.arange(3,11),styles[1], markersize=10, linewidth = 3.3)
  p3=plt.plot(a[:,0],a[:,3]/2**np.arange(3,11),styles[2], markersize=10, linewidth = 3.3)
  plt.ylim((-1.5, -0.95))
  #p6=plt.plot(a8[:,0]/2**8,a8[:,1]/2**8,styles[5], markersize=8, linewidth = 3.3)
  #p7=plt.plot(a9[:,0]/2**9,a9[:,1]/2**9,styles[6], markersize=8, linewidth = 3.3)
  plt.grid('on')
  plt.ylabel('Lowest Energy Density', fontsize = 16)
  plt.xlabel(' log(N)', fontsize = 16)
  plt.title(titlename, fontsize = 16)
  plt.legend([p1[0],p2[0],p3[0]],('HN3', 'HN5', 'HNNP'),'center')
  plt.savefig(filename)


  
  return 0;
  
  


if __name__ == '__main__':
  main()

