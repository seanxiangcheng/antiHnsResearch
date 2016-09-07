import sys
import numpy as np
import random
import matplotlib.pyplot as plt


def main():
  lab = np.loadtxt("Results_lab0z.txt")
  lab_sr = np.loadtxt("Results_sr_lab0z.txt")
  cheetah = np.loadtxt("Results_cheetah.txt")
  cheetah_sr = np.loadtxt("Results_sr_cheetah.txt")
  styles = ['-bh', '-gv', '-rd', '-cs', '-ko']
  styles2 = ['--bh', '--gv', '--rd', '--cs', '--ko']
  #plt.figure(1)
  plt.figure(1, figsize=(8,8))
  for i in range(lab.shape[1]-1):
    plt.loglog(lab[:,0], lab[:,i+1], styles[i], markersize=10, linewidth = 3.3)
  
  for i in range(lab.shape[1]-1):
    plt.loglog(lab_sr[:,0], lab_sr[:,i+1], styles2[i], markersize=10, linewidth = 3.3)
  plt.title('Lab: MPI_Bcast (solid line) vs Handcoded (dashed line)', fontsize = 16)
  plt.grid(True)  
  plt.xlim([1.84, 136.8])
  plt.ylim([np.amin(lab[:,1:])*0.80, np.amax(lab_sr[:,1:])*1.29])
  plt.xlabel('Number of processes', fontsize = 16);
  plt.ylabel('Time spent (s)', fontsize = 16)
  plt.xticks([2, 4, 8, 16, 32, 64, 128], [2, 4, 8, 16, 32, 64, 128], fontsize = 16)
  plt.yticks(fontsize = 16)
  plt.savefig('lab_plot.pdf')




  plt.figure(2, figsize=(8,8))
  for i in range(cheetah.shape[1]-1):
    plt.loglog(cheetah[:,0], cheetah[:,i+1], styles[i], markersize=10, linewidth = 3.3)
  
  for i in range(cheetah_sr.shape[1]-1):
    plt.loglog(cheetah_sr[:,0], cheetah_sr[:,i+1], styles2[i], markersize=10, linewidth = 3.3)
  plt.title('Cheetah: MPI_Bcast (solid line) vs Handcoded (dashed line)', fontsize = 16)
  plt.grid(True)  
  plt.xlim([1.84, 136.8])
  plt.ylim([np.amin(cheetah[:,1:])*0.80, np.amax(cheetah_sr[:,1:])*1.29])
  plt.xlabel('Number of processes', fontsize = 16)
  plt.ylabel('Time spent (s)', fontsize = 16)
  plt.xticks([2, 4, 8, 16, 32, 64, 128], [2, 4, 8, 16, 32, 64,128], fontsize = 16)
  plt.yticks(fontsize = 16)
  plt.savefig('cheetah_plot.pdf')
  
  #plt.legend([plot[0], plot[1], plot[2], plot[3], plot[4]],('1 Byte', '1K', '4K', '64K', '1M'))
    #plt.set_markersize(3)

  #plt.show()
  
  return 0;
  
  


if __name__ == '__main__':
  main()
