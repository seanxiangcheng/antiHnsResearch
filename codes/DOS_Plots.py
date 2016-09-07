import sys
import numpy as np
import matplotlib.pyplot as plt

def main():
  titlename = 'HN5 Anti-Ferr Ising'
  filename = 'H5DOS.pdf'
  a3 = np.loadtxt("H5AFDOS3_r0")
  a4 = np.loadtxt("H5AFDOS4_r0")
  a5 = np.loadtxt("H5AFDOS5_r0")
  a6 = np.loadtxt("H5AFDOS6_r0")
  a7 = np.loadtxt("H5AFDOS7_r0")
  a8 = np.loadtxt("H5AFDOS8_r0")
  a9 = np.loadtxt("H5AFDOS9_r0")
  styles = ['-bh', '-mo', '-gd','-y', '-r', '-c', '-k']
  
  print "%-5d %-6d" % (2**3, a3[-1,0]) 
  print "%-5d %-6d" % (2**4, a4[-1,0]) 
  print "%-5d %-6d" % (2**5, a5[-1,0]) 
  print "%-5d %-6d" % (2**6, a6[-1,0]) 
  print "%-5d %-6d" % (2**7, a7[-1,0]) 
  print "%-5d %-6d" % (2**8, a8[-1,0]) 
  print "%-5d %-6d" % (2**9, a9[-1,0]) 

  plt.figure(1, figsize=(8,8))
  p1=plt.plot(a3[:,0]/2**3,a3[:,1]/2**3,styles[0], markersize=8, linewidth = 3.3)
  p2=plt.plot(a4[:,0]/2**4,a4[:,1]/2**4,styles[1], markersize=8, linewidth = 3.3)
  p3=plt.plot(a5[:,0]/2**5,a5[:,1]/2**5,styles[2], markersize=8, linewidth = 3.3)
  p4=plt.plot(a6[:,0]/2**6,a6[:,1]/2**6,styles[3], markersize=8, linewidth = 3.3)
  p5=plt.plot(a7[:,0]/2**7,a7[:,1]/2**7,styles[4], markersize=8, linewidth = 3.3)
  p6=plt.plot(a8[:,0]/2**8,a8[:,1]/2**8,styles[5], markersize=8, linewidth = 3.3)
  p7=plt.plot(a9[:,0]/2**9,a9[:,1]/2**9,styles[6], markersize=8, linewidth = 3.3)
  plt.grid('on')
  plt.ylim((-0.01,np.amax(a9[:,1]/2**9)*1.02))
  plt.xlim((np.amin(a9[:,0]/2**9))*1.018)
  plt.xlabel('Energy Density', fontsize = 16)
  plt.ylabel('Log(g) / N', fontsize = 16)
  plt.title(titlename, fontsize = 16)
  plt.legend([p1[0],p2[0],p3[0],p4[0],p5[0],p6[0], p7[0]],('N=8', 'N=16', 'N=32', 'N=64', 'N=128', 'N=256', 'N=512'),'upper right')
  #plt.savefig(filename)
  plt.show()


  
  return 0;
  
  


if __name__ == '__main__':
  main()

