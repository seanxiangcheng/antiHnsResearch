# Python code: HN3 and HN5 Phase diagram
import numpy as np
import matplotlib.pyplot as plt
import HNRG as hnrg

def main():
  nt = 100
  ms = np.linspace(0.001, 0.5, num = nt)
  ts = np.array([1.0, 1.5])
  n = 1000
  l = 500
  bins = np.linspace(-1.5, 1.5, num=100);
  ys = np.array([-2., -1.5, -1., -.5, 0, 0.5, 1., 1.5, 2.])
  y = 0
  styles = ['b-o', 'k-x', 'c->', 'r-s', 'g-d']
  plt.figure(1, figsize = (10,18))
  hp = hnrg.hp6()
  #ks = np.zeros((num_init,))
  for yi in range(ys.size):
    y = ys[yi]
    plt.subplot(int(np.ceil(ys.size/2.)), 2, yi+1)
    for i in range(ts.size):
      t = ts[i]
      m = np.exp(2.0/t)
      js = -np.log(hp.rg_n_k(m, y, n, l))*(1.0)*t/4.0
      hist, centers = np.histogram(js, bins = bins)
      centers = (centers[:-1]+centers[1:])/2.
      hist = hist/1.0/np.sum(hist)
      plt.plot(centers, hist, styles[i], label = 'T='+str(ts[i]))
  #plt.semilogy(np.linspace(0.3333, 1, num = nt), np.ones((nt,)), 'k-', linewidth = 2)
    plt.grid('on')
    plt.title('y='+str(y))
    if yi==0 or yi==4:
      plt.legend(loc='upper center')
    elif yi==1 or yi==2 or yi==3 or yi==5 or yi==6:
      plt.legend(loc = 'upper left')
    else:
      plt.legend(loc = 'upper right')


    #plt.xlim([0., 1.01])
    #plt.ylim([0, 1.01])
    #plt.legend(loc = 2)
    #plt.ylabel('$J$', fontsize = 16)
    if yi==ys.size-1:
      plt.xlabel('$J$', fontsize = 16)
    #plt.title('HP6,y='+str(y)+',steps:'+str(n-l)+'~'+str(n), fontsize = 16)
  #plt.savefig('HP6_AFM_J_dist_np.pdf')
  #plt.savefig('HP6_AFM_J_dist_np.png', dpi = 100)
  
  jfps = np.zeros(())
  for yi in range(ys.size):
    y = ys[i]
    plt.subplot(int(np.ceil(ys.size/2.)), 2, yi+1)
    for i in range(ts.size):
      t = ts[i]
      m = np.exp(2.0/t)
      js = -np.log(hp.rg_n_k(m, y, n, l))*(1.0)*t/4.0
      hist, centers = np.histogram(js, bins = bins)
      centers = (centers[:-1]+centers[1:])/2.
      hist = hist/1.0/np.sum(hist)
      plt.plot(centers, hist, styles[i], label = 'T='+str(ts[i]))
  plt.show()
  

    
    

if __name__ == '__main__':
  main()
  
