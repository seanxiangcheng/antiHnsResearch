# Python code: HN3 and HN5 Phase diagram
import numpy as np
import matplotlib.pyplot as plt

def main():
  N = 50
  lw = 2.5
  fs = 16
  ms = 4
  styles = ['-bh', '-c>', '-y^', '-rd','-m>', '-b*', '-gs','-r8','-ko']
  mu = np.linspace(0.001, 1, num = N)
  kappa = np.zeros((N,))
  y = np.array([1.5, 1, 0.5, 0, -0.5, -1, -1.5, -2, -10])
  
  plt.figure(figsize = (10,10))
  for i in range(y.size):
    for j in range(mu.size):
      kappa[j] = kap_mu(1.0/mu[j], y[i])
    plt.plot(mu, kappa, styles[i], linewidth = lw, markersize = ms, label = 'y='+str(y[i]))
  plt.plot(mu, np.zeros((N,)), '--k', linewidth = 3, label = '$\kappa = 0$')
  plt.grid('on')
  plt.xlim([0., 1.001])
  plt.ylim([-0.45, 2.23])
  plt.ylabel('$\kappa$', fontsize = fs*1.5)
  plt.xlabel('$1/\mu$ ', fontsize = fs*1.2)
  plt.title('y='+str(y), fontsize = fs)
  plt.title('HN networks with different y', fontsize = fs)
  plt.legend(loc = "upper left")
  plt.savefig('AFM_HN_ys_kvsmu.pdf')
  #plt.show()


# equation to find fixed point kappa; 
def kap_mu(m, y):
  value = 0.50*(m**(1+y) + m**y -2*m + np.sqrt(4.0*(m**(1.0+y)+m**y-1.0)+(m**(1.0+y)+m**y-2.0*m)**2))
  return value



if __name__ == '__main__':
  main()
  
