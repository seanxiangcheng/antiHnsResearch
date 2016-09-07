# Python code: HN3 and HN5 Phase diagram
import numpy as np
import matplotlib.pyplot as plt

N = 250   # number of RG iterations

# recursive equation to find lambda; (from recursive equation)
def lam_rec(m, k, y):
  value = m**(2.0*y)*(1+m)*(1+k)**2/(2*(1+2*m*k+k**2))
  return value


# recursive equation to find kappa; (from recursive equation)
def kap_rec(m, k, lam):
  value = 2*k*lam*(1+m)/(1+2*m*k+k**2)
  return value


# equation to find fixed point kappa; 
#(from fixed point solutions bsaed on the recursive equations)
# DONOT USE it when y != 1
def kap_mu(m):
  value = 0.50 * (-m + m**2 + np.sqrt(-4+4*m+5*m**2-2*m**3+m**4))
  return value


# equation to find fixed point lambda; 
# (from fixed point solutions bsaed on the recursive equations)
# DONOT USE it when y != 1
def lam_mu(m):
  value = m/4.0 * (-m + 2 + m**2 + np.sqrt(-4+4*m+5*m**2-2*m**3+m**4))
  return value
  

def main():
  C = np.zeros((N,))
  kappa = np.zeros((N,))
  lambda_rg = np.zeros((N,))
  num_init = 24
  kl_end = np.zeros((num_init, 2))
  kappa_end = np.zeros((num_init,))
  lambda_end = np.zeros((num_init,))
  mu = np.linspace(0.01, 1, num = num_init)
  y = 1     # HN3: y=0; HN5: y=1
  C[0] = 1.0
  
  #plt.figure(1, figsize = (6,12))
  plt.figure(figsize = (12,6))
  plt.subplot(1,2, 1)
  for k in range(mu.size):
    m = 1/mu[k]
    kappa[0] = m**2
    lambda_rg[0] = m**(2*y)
    #kappa_end[k] = 1/kap_mu(m)
    #lambda_end[k] = 1/lam_mu(m) # only use it when y = 1 
    for i in range(N)[1:]:
      lambda_rg[i] = lam_rec(m, kappa[i-1], y)
      kappa[i] = kap_rec(m, kappa[i-1], lambda_rg[i-1])
    #kappa = 1/kappa
    #lambda_rg = 1/lambda_rg
    kl_end[k,0] = lambda_rg[-1]
    kl_end[k,1] = kappa[-1]
    plt.semilogy(lambda_rg, kappa,'-ok', markersize = 5)
    plt.semilogy(lambda_rg[0], kappa[0],'-<b', markersize = 9)
    plt.semilogy(lambda_rg[-1], kappa[-1],'-or', markersize = 9)
  plt.semilogy(kl_end[:,0], kl_end[:,1], 'r', linewidth = 2.5)
  print "kappa[ %.4e, %.4e ], lambda [%.4e, %.4e]" % (np.min(kl_end[:,1]), np.max(kl_end[:,1]), np.min(kl_end[:,0]), np.max(kl_end[:,0]))
  #plt.plot(lambda_end, kappa_end, '-*g', markersize = 10, linewidth = 2) # plot using fixed point solutions of both kappa and lambda
  plt.grid('on')
  #plt.xlim([0., 1.01])
  #plt.ylim([0, 1.01])
  plt.xlabel('$1/\lambda$', fontsize = 18)
  plt.ylabel('$\kappa$ ', fontsize = 18)
  plt.title('y='+str(y), fontsize = 18)
  #plt.show()
  
  plt.subplot(1,2, 2)
  plt.plot(mu, kl_end[:,1],'-ro', markersize = 8, linewidth = 2.5, label = ('y='+str(y)+' at 1K RG steps'))
  #plt.plot(mu, kappa_end,'-bs', markersize = 8, linewidth = 2.5, label = 'y=1 using fixed point sol')
  plt.grid('on')
  #plt.xlim([0., 1.01])
  #plt.ylim([0, 1.01])
  #plt.legend(loc = "upper right")
  plt.ylabel('$\kappa$', fontsize = 18)
  plt.xlabel('$1/\mu$', fontsize = 18)
  plt.title('y='+str(y), fontsize = 18)
  #plt.savefig('AFM_HN_yn'+str(-y)+'.pdf')
  plt.show()
  
if __name__ == '__main__':
  main()
  
