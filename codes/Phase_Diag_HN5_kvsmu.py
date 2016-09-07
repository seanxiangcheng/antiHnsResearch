# Python code: HN3 and HN5 Phase diagram
import numpy as np
import matplotlib.pyplot as plt

N = 1000   # number of RG iterations

def lam_rec(m, k, y):
  value = m**(2*y)*(1+m)*(1+k)**2/(2*(1+2*m*k+k**2))
  return value
  
def kap_rec(m, k, lam):
  value = 2*k*lam*(1+m)/(1+2*m*k+k**2)
  return value

def kap_mu(m):
  value = 0.50 * (-m + m**2 + np.sqrt(-4+4*m+5*m**2-2*m**3+m**4))
  return value

def main():
  num_init = 20
  kappa = np.zeros((N,))
  kappa_mu = np.zeros((num_init,))
  lambda_rg = np.zeros((N,))
  kappa_end = np.zeros((num_init, 1))
  mu = np.linspace(0.5*(np.sqrt(5.0)-1), 1, num = num_init)
  print mu
  y = 1     # HN3: y=0; HN5: y=1
  
  plt.figure(1, figsize = (6,6))
  for k in range(mu.size):
    m = mu[k]
    kappa[0] = m**2
    lambda_rg[0] = m**(2*y)
    kappa_mu[k] = kap_mu(m)
    for i in range(N)[1:]:
      lambda_rg[i] = lam_rec(m, kappa[i-1], y)
      kappa[i] = kap_rec(m, kappa[i-1], lambda_rg[i-1])
    #kappa = 1/kappa
    #kappa_end[k] = kappa[-1]
  #kappa_mu = 1/kappa_mu
  plt.plot(mu, kappa_end,'or', markersize = 6)
  plt.plot(mu, kappa_mu,'b', linewidth = 3)
  plt.grid('on')
  #plt.xlim([0., 1.01])
  #plt.ylim([0, 1.01])
  plt.xlabel('$\mu$ or $1/\mu$', fontsize = 18)
  plt.ylabel('$\kappa$ or $1/\kappa$', fontsize = 18)
  plt.title('HN5, AFM, 1000 steps, Phase diagram', fontsize = 18)
  #plt.savefig('AFM_HN5_kappavsmu.pdf')
  plt.show()
  
if __name__ == '__main__':
  main()
  
