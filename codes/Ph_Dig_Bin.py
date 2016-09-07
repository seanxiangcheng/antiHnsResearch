# Python code: HN3 and HN5 Phase diagram
import numpy as np
import matplotlib.pyplot as plt

Nmin = 1000
N = 5000   # number of RG iterations

# recursive equation to find lambda; (from recursive equation)
def lam_rec(m, k):
  value = ((m+k)**2.0)/((1+m*k)**2.0)
  return value


# recursive equation to find kappa; (from recursive equation)
def kap_rec(m, k, lam):
  value = ((1+m)**2)*k*lam/((1+m*k)**2)
  return value



def main():
  kn = 0.0
  ko = 0.0
  ln = 0.0
  lo = 0.0
  num_init = 50
  ks = np.zeros((N - Nmin,))
  mu = 0.3
  
  plt.figure(1, figsize = (8,6))
  m =1/ mu
  ko = m**2
  lo = 1
  for i in range(N):
    ln = lam_rec(m, ko)
    kn = kap_rec(m, ko, lo)
    lo = ln
    ko = kn
    if i >= Nmin:
      ks[i-Nmin] = kn
  hist, bins = np.histogram(ks, bins = 100)
  centers = (bins[:-1] + bins[1:])/2
  print 'hist\n', hist
  print 'bins\n', bins
  plt.plot(centers, hist, '-o')
  plt.grid('on')
  #plt.xlim([0., 1.01])
  #plt.ylim([0, 1.01])
  #plt.legend(loc = 2)
  plt.ylabel('visits', fontsize = 18)
  plt.xlabel('$\kappa$', fontsize = 18)
  plt.title("mu="+str(mu)+'; N='+str(N)+';Nmin='+str(Nmin), fontsize = 18)
  #plt.savefig('Histo_1.pdf')
  plt.show()
  
if __name__ == '__main__':
  main()
  
