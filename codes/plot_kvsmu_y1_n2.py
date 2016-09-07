# Python code: HN3 and HN5 Phase diagram
import numpy as np
import matplotlib.pyplot as plt

N = 1000   # number of RG iterations

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
  num_init = 20
  kl_end = np.zeros((num_init, 2))
  kappa_end = np.zeros((num_init,))
  lambda_end = np.zeros((num_init,))
  mu = np.linspace(0.02, 1, num = num_init)
  print mu
  styles = ['-ro','-g>', '-bd','-ks','-go','-b>', '-kd','-rs']
  y = [1.5, 1, 0.5, 0, -0.5, -1, -1.5, -2.0]     # HN3: y=0; HN5: y=1
  C[0] = 1.0
  
  plt.figure(figsize = (6,6))
  for yi in range(len(y)):
    for k in range(mu.size):
      m = 1/mu[k]
      kappa[0] = m**2
      lambda_rg[0] = m**(2*y[yi])
      for i in range(N)[1:]:
        lambda_rg[i] = lam_rec(m, kappa[i-1], y[yi])
        kappa[i] = kap_rec(m, kappa[i-1], lambda_rg[i-1])
    #kappa = 1/kappa
    #lambda_rg = 1/lambda_rg
      kl_end[k,0] = lambda_rg[-1]
      kl_end[k,1] = kappa[-1]
    plt.plot(mu, kl_end[:,1], styles[yi])
  plt.grid('on')
  plt.xlim([0., 1.01])
  plt.ylim([0, 2.01])
  #plt.legend(loc = 2)
  plt.ylabel('$\kappa$ ', fontsize = 18)
  plt.xlabel('$1/\mu$', fontsize = 18)
  plt.title('HN5, AFM, y=1', fontsize = 18)
  #plt.savefig('AFM_HN_y8_kappavsmu.pdf')
  plt.show()
  
if __name__ == '__main__':
  main()
  
