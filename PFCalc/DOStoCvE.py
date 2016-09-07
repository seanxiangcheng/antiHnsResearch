import numpy as np
import matplotlib.pyplot as plt

ACCU = 20


def main():
  ### variables for consistence ####
  fn_ds = "H3AFDOS6_r0"
  nt = 100
  t_min = 0.1
  t_max = 10
  style = '-'
  ms = 3            # markersize
  lw = 2            # linewidth
  fs = 12           # fontsize
  ts = np.linspace(t_min, t_max, num = nt) 
  energy = np.zeros((ts.size))  # internal energy
  cv = np.zeros((ts.size))      # specific heat; it is <E^2> at first
  s = np.zeros((ts.size))       # entropy; it is P*Z at first
  fe = np.zeros((ts.size))      # free energy
  
  ### load data ###
  dos = np.loadtxt(fn_ds)
  #es = dos[:,0]
  #logg = dos[:,1]
  
  ### Calculations ###
  for i in range(ts.size):
    t = ts[i]
    dos_t = DOS_t(dos, t)
    es = dos_t[:,0]
    logg = dos_t[:,1]
    z = 0.0
    for j in range(es.size):
      temp_exp = np.exp(logg[j]-es[j]/t)
      z = z + temp_exp
      energy[i] = energy[i] + es[j]*temp_exp
      cv[i] = cv[i] + es[j]*es[j]*temp_exp

    fe[i] = -t*np.log(z)
    energy[i] = energy[i]/z
    cv[i] = (cv[i]/z - energy[i]*energy[i])/t/t
    for j in range(es.size):
      temp_exp = np.exp(logg[j]-es[j]/t)
      s[i] = s[i] - temp_exp/z * np.log(temp_exp/z)
  #print 'energy\n', energy
  #print 'cv\n', cv
  #print 'free e\n', fe
  #print 's\n', s
  
  plt.figure(1, figsize = (8,6))
  
  plt.subplot(221)
  plt.plot(ts, energy, style, markersize = ms, linewidth = lw)
  plt.xlabel('Temperature', fontsize = fs)
  plt.ylabel('Internal Energy', fontsize = fs)
  
  plt.subplot(222)
  plt.plot(ts, fe, style, markersize = ms, linewidth = lw)
  plt.xlabel('Temperature', fontsize = fs)
  plt.ylabel('Free Energy', fontsize = fs)
  
  plt.subplot(223)
  plt.plot(ts, cv, style, markersize = ms, linewidth = lw)
  plt.xlabel('Temperature', fontsize = fs)
  plt.ylabel('Specific Heat', fontsize = fs)
  
  plt.subplot(224)
  plt.plot(ts, s, style, markersize = ms, linewidth = lw)
  plt.xlabel('Temperature', fontsize = fs)
  plt.ylabel('Entropy', fontsize = fs)
  
  
  
  plt.show()

    




def DOS_t(dos, t):
  temp = 0.0
  max_ex = 0.0
  for i in range(dos.shape[0]):
    temp = dos[i,1]-dos[i, 0]/t
    if temp > max_ex:
      max_ex = temp
  new_dos = np.zeros(dos.shape)
  count = 0
  for i in range(dos.shape[0]):
    temp = dos[i,1]-dos[i, 0]/t
    if temp>(max_ex - ACCU):
      new_dos[count,:] = dos[i,:]
      count = count+1
  return new_dos[:count, :]



      
  



















if __name__ == '__main__':
  main()
