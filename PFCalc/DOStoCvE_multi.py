import numpy as np
import matplotlib.pyplot as plt

ACCU = 25
Styles = ['-b', '-r', '-g','-k']
Filenames = ['HPAFDOS6_r0', 'HPAFDOS7_r0', 'HPAFDOS8_r0', 'HPAFDOS9_r0']
Ns = np.array([64, 128,256,512])
def main():
  plt.figure(1, figsize = (12,12))
  for i in range(len(Filenames)):
    Calc_plot(Filenames[i], i)
  plt.subplot(221)
  plt.legend(loc = "upper left")
  
  data = np.loadtxt("HPAF_E_109")
  plt.subplot(221)
  plt.plot(data[:,0], data[:,1:], '--')
  plt.xlim([0, 1.2])
  plt.ylim([-1.5, -1.25])
  plt.savefig(Filenames[1][:4] + 'e_fe_cv_s.pdf')
  plt.show()

def Calc_plot(filename, sty_i):
  ### variables for consistence ####
  fn_ds = filename
  style = Styles[sty_i%len(Styles)]
  N = Ns[sty_i]
  nt = 100
  t_min = .01
  t_max = 1.2
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
    threshd = logg[0] - es[0]/t
    for j in range(es.size):
      temp_exp = np.exp(logg[j]-es[j]/t - threshd )
      z = z + temp_exp
      energy[i] = energy[i] + es[j]*temp_exp
      cv[i] = cv[i] + es[j]*es[j]*temp_exp

    fe[i] = -t*(np.log(z)+threshd)/N
    energy[i] = energy[i]/z/N
    cv[i] = (cv[i]/z - energy[i]*N*energy[i]*N)/t/t/N
    #if N==512 and t < 1:
      #print t, ":", energy[i]
    for j in range(es.size):
      temp_exp = np.exp(logg[j]-es[j]/t - threshd)
      s[i] = s[i] - temp_exp/z * np.log(temp_exp/z)
  #s = s/N
  #print 'energy\n', energy
  #print 'cv\n', cv
  #print 'free e\n', fe
  #print 's\n', s
    
  plt.subplot(221)
  plt.plot(ts, energy, style, markersize = ms, linewidth = lw, label='N='+str(N))
  plt.xlabel('Temperature', fontsize = fs)
  plt.ylabel('Internal Energy Density', fontsize = fs)
  
  plt.subplot(222)
  plt.plot(ts, fe, style, markersize = ms, linewidth = lw)
  plt.xlabel('Temperature', fontsize = fs)
  plt.ylabel('Free Energy Density', fontsize = fs)
  
  plt.subplot(223)
  plt.plot(ts, cv, style, markersize = ms, linewidth = lw)
  plt.xlabel('Temperature', fontsize = fs)
  plt.ylabel('Specific Heat Density', fontsize = fs)
  #plt.xlim([0.35, 0.62])
  #plt.ylim([0.12, 0.34])
  
  plt.subplot(224)
  plt.plot(ts, s/N, style, markersize = ms, linewidth = lw)
  plt.xlabel('Temperature', fontsize = fs)
  plt.ylabel('Entropy Density', fontsize = fs)
  
  
  

    




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
