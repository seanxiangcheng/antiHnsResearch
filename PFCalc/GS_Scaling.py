import numpy as np
import matplotlib.pyplot as plt




def main():
  styles = ['-bo', '-rs', '-gh','-c>','-kv']
  filenames = ['H3AF_E_109', 'H3AF_E_110', 'H3AF_E_112','H3AF_E_116', 'H3AF_E_117']
  Ns = 2**np.array([9, 10, 12, 16, 17])
  dt = 0.001/2**np.arange(10)
  lw = 2.5
  fs = 15
  plt.figure(1, figsize = (8,6))
  for i in range(len(filenames)):
    data = np.loadtxt(filenames[i])
    plt.loglog(dt, 1 + data[-1,1:], styles[i], linewidth = 2.5, label = 'N='+str(Ns[i]))
  plt.legend(loc = 'lower center')
  plt.xlabel('dT', fontsize = fs)
  plt.ylabel('$E(T=0.01) - (-1)$', fontsize = fs)
  plt.title('HN3 lowest energy at T=0.01 with differnt dT', fontsize = fs)
  plt.ylim([0.004, 0.04])
  plt.xlim([np.amin(dt)*0.8, np.amax(dt)*1.2])
  plt.yticks([0.004, 0.006, 0.008, 0.01, 0.02, 0.03, 0.04], [0.004, 0.006, 0.008, 0.01, 0.02, 0.03, 0.04])

  plt.savefig(filenames[0][:2]+'lwe_dts.pdf')
  plt.show()




if __name__ == '__main__':
  main()
  
