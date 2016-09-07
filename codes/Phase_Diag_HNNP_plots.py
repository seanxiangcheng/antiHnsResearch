# Python code: HN3 and HN5 Phase diagram
import numpy as np
import matplotlib.pyplot as plt

N = 1000   # number of RG iterations

def lam_rec(m, k, y):
  value = m**2*(1+m)*(1+k)**2/(2*(1+2*m*k+k**2))
  return value
  
def kap_rec(m, k, lam):
  value = 2*(1+m)*k*lam/(1+2*m*k+k**2)
  return value

def kap_mu(m):
  value = 0.50 * (-m + m**2 + np.sqrt(-4+4*m+5*m**2-2*m**3+m**4))
  return value

def main():
  C = np.zeros((N,))
  kappa = np.zeros((N,))
  lambda_rg = np.zeros((N,))
  num_init = 20
  kl_end = np.zeros((num_init, 2))
  kappa_end = np.zeros((num_init,))
  mu = np.linspace(0.02, 1, num = num_init)
  print mu
  y = 1     # HN3: y=0; HN5: y=1
  C[0] = 1.0
  
  plt.figure(1, figsize = (6,6))
  for k in range(mu.size):
    m = 1/mu[k]
    kappa[0] = m**2
    lambda_rg[0] = m**2
    kappa_end[k] = 1/kap_mu(m)
    for i in range(N)[1:]:
      lambda_rg[i] = lam_rec(m, kappa[i-1], y)
      kappa[i] = kap_rec(m, kappa[i-1], lambda_rg[i-1])
    kappa = 1/kappa
    lambda_rg = 1/lambda_rg
    kl_end[k,0] = lambda_rg[-1]
    kl_end[k,1] = kappa[-1]
    plt.plot(lambda_rg, kappa,'-ok', markersize = 5)
    plt.plot(lambda_rg[0], kappa[0],'-<b', markersize = 9)
    plt.plot(lambda_rg[-1], kappa[-1],'-or', markersize = 9)
  plt.plot(kl_end[:,0], kappa_end,'b', linewidth = 2)
  plt.grid('on')
  plt.xlim([0., 1.01])
  plt.ylim([0, 1.01])
  plt.xlabel('$\lambda$ or $1/\lambda$', fontsize = 18)
  plt.ylabel('$\kappa$ or $1/\kappa$', fontsize = 18)
  plt.title('HN5, AFM, 1000 steps, blue line<-solution', fontsize = 18)
  #plt.savefig('AFM_HNNP_kappavslambda.pdf')
  #plt.show()
  
  plt.figure(2, figsize = (6,6))
  plt.plot(mu, kl_end[:,1],'bo', markersize = 8)
  plt.plot(mu, kappa_end,'-r', linewidth = 3)
  plt.grid('on')
  plt.xlim([0., 1.01])
  plt.ylim([0, 1.01])
  plt.ylabel('$\kappa$ or $1/\kappa$', fontsize = 18)
  plt.xlabel('$\mu$ or $1/\mu$', fontsize = 18)
  plt.title('HN5, AFM, 1000 steps', fontsize = 18)
  #plt.savefig('AFM_HNNP_kappavsmu.pdf')
  plt.show()
  
if __name__ == '__main__':
  main()
  
