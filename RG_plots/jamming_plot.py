import numpy as np
import matplotlib.pyplot as plt

files = ["HPAF_E_109", "HPAF_E_110", "HPAF_E_112", "HPAF_E_114"]
gs = [-1514., -6062., -24252.]
Ls = [2**10, 2**12, 2**14]
dts = 0.001/np.power(2, range(10))

plt.figure(1, figsize = (10,12))
for i in range(len(gs)):
  plt.subplot(2,2, i+1)
  data = np.loadtxt(files[i])
  plt.plot(data[:,0], (data[:,1:]-gs[i]/Ls[i]))
  #print 'sa:', data[-1,:]
  #print 'gs:', gs[i]/Ls[i]
  plt.xlim([0,1])
  plt.ylim([0, 0.15])
  plt.xlabel("T")
  plt.ylabel("E - GS")
  plt.title('L='+str(Ls[i]))
  
plt.subplot(2,2,4)
plt.loglog(dts, (data[-1,1:]-gs[i]/Ls[i]), 'o-')
print dts
print (data[-1,1:]-gs[i]/Ls[i])
#z = np.linalg.lstsq(np.array([dts, np.ones(10)]), (data[-1,1:]-gs[i]/Ls[i]))[0]
#plt.loglog(dts, z[0]*dts+z[1], '-')
plt.xlim([np.amin(dts)*0.9, np.amax(dts)*1.1])
plt.ylim([0.0036, 0.022])
#plt.savefig("HP_jamming_plot.pdf")
plt.show()
