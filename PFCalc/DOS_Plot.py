import numpy as np
import matplotlib.pyplot as plt




def main():
  styles = ['-b', '-r', '-g','-k']
  filenames = ['H3AFDOS9_r0', 'H5AFDOS9_r0', 'HPAFDOS9_r0']
  N = 512
  lw = 2.5
  fs = 15
  plt.figure(1, figsize = (8,6))
  for i in range(len(filenames)):
    data = np.loadtxt(filenames[i])
    plt.plot(data[:,0]/N, data[:,1]/N, styles[i], linewidth = 2.5, label = filenames[i][:4])
  plt.legend(loc = 'lower center')
  plt.xlabel('Energy density', fontsize = fs)
  plt.ylabel('$\log[DOS]/L$', fontsize = fs)
  plt.title('HNs N = 512, Density of states from Wang-Landau', fontsize = fs)
  plt.savefig('DOSs.pdf')
  plt.show()




if __name__ == '__main__':
  main()
  
