# Python code: HN3 and HN5 Phase diagram
import numpy as np
import matplotlib.pyplot as plt

N = 1000   # number of RG iterations

# recursive equation to find lambda; (from recursive equation)
def lam_rec(m, k):
  value = ((m+k)**2.0)/((1+m*k)**2.0)
  return value


# recursive equation to find kappa; (from recursive equation)
def kap_rec(m, k, lam):
  value = ((1+m)**2)*k*lam/((1+m*k)**2)
  return value



def main():
  C = np.zeros((N,))
  kappa_rg = np.zeros((N,))
  lambda_rg = np.zeros((N,))
  num_init = 50
  kl_end = np.zeros((num_init, 2))
  mu = np.linspace(0.01, 1, num = num_init)
  C[0] = 1.0
  
  plt.figure(1, figsize = (8,6))
  #plt.subplot(121)
  for k in range(mu.size):
    m =1/ mu[k]
    kappa_rg[0] = m**2
    lambda_rg[0] = 1
    for i in range(N)[1:]:
      lambda_rg[i] = lam_rec(m, kappa_rg[i-1])
      kappa_rg[i] = kap_rec(m, kappa_rg[i-1], lambda_rg[i-1])
    kappa_rg = kappa_rg
    lambda_rg = lambda_rg
    kl_end[k,0] = lambda_rg[-1]
    kl_end[k,1] = kappa_rg[-1]
    #plt.loglog(lambda_rg, kappa_rg,'-ok', markersize = 5)
    #plt.loglog(lambda_rg[0], kappa_rg[0],'<b', markersize = 9)
    #plt.loglog(lambda_rg[-1], kappa_rg[-1],'or', markersize = 7)
  #plt.plot(kl_end[:,0], kl_end[:,1], 'r', linewidth = 2, label = 'Fixed point solution')
  #plt.plot(lambda_end, kappa_end, '-*g', markersize = 10, linewidth = 2) # plot using fixed point solutions of both kappa and lambda
  #plt.grid('on')
  #plt.xlim([0., 1.01])
  #plt.ylim([0, 1.01])
  #plt.legend(loc = 2)
  #plt.xlabel('$\lambda$', fontsize = 18)
  #plt.ylabel('$\kappa$', fontsize = 18)
  #plt.title('HNNP, AFM, $\lambda^{(0)}=1, \kappa^{(0)} = \mu^2$'+'RG'+str(N), fontsize = 18)
  #plt.savefig('AFM_HNNP_kappavslambda.pdf')
  
  #plt.subplot(122)
  plt.semilogy(mu, kl_end[:,1],'-ro', markersize = 5)
  plt.grid('on')
  #plt.xlim([0., 1.01])
  #plt.ylim([0, 1.01])
  #plt.legend(loc = 2)
  plt.ylabel('$\kappa$', fontsize = 18)
  plt.xlabel('$1/\mu$', fontsize = 18)
  plt.title('HNNP, AFM, $\lambda^{(0)}=1, \kappa^{(0)} = \mu^2$, RG step'+str(N), fontsize = 18)
  #plt.savefig('AFM_HNNP_dig'+str(N)+'_kvsmu.pdf')
  plt.show()
  
if __name__ == '__main__':
  main()
  
