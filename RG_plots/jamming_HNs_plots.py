import numpy as np
import matplotlib.pyplot as plt

files = ["H3AF_E_114", "H5AF_E_114", "HPAF_E_114", "H6AF1D_E_k14_r0"]
Efiles = ["H3AF_E_low114", "H5AF_E_low114", "HPAF_E_low114", "H6AF1D_E_low_k14_r0"]

gs = [-16384., -19112., -24252., -21064.]
ylims = [[0, 0.14], [0,0.1], [0, 0.1], [0, 0.1]]
Ls = [2**14, 2**14, 2**14, (2**14+1)]
t_loc = [[0.05, 0.123], [0.05, 0.0883], [0.05, 0.0883], [0.05, 0.0883]]
strs = ["HN3", "HN5", "HNNP", "HN6"]
stls = ['b-o', 'g->', 'r-s', 'k-d']
dts = 0.001/np.power(2, range(10))

plt.figure(1, figsize = (20,5.1))
for i in range(len(gs)):
  plt.subplot(1,4, i+1)
  data = np.loadtxt(files[i])
  plt.plot(data[:,0], (data[:,1:]-gs[i]/Ls[i]))
  #print 'sa:', data[-1,:]
  #print 'gs:', gs[i]/Ls[i]
  plt.xlim([0,0.8])
  plt.ylim(ylims[i])
  plt.xlabel("$T$", fontsize = 16)
  if i==0:
    plt.ylabel("$E(T) - GS$", fontsize = 16)
  plt.text(t_loc[i][0], t_loc[i][1], strs[i], fontsize = 20)
  plt.text(t_loc[i][0], t_loc[i][1]*0.9, 'N=16384', fontsize = 20)
  plt.grid('on')
  #plt.title('L='+str(Ls[i]))
  
#plt.subplot(2,2,4)
#plt.loglog(dts, (data[-1,1:]-gs[i]/Ls[i]), 'o-')
#print dts
#print (data[-1,1:]-gs[i]/Ls[i])
#z = np.linalg.lstsq(np.array([dts, np.ones(10)]), (data[-1,1:]-gs[i]/Ls[i]))[0]
#plt.loglog(dts, z[0]*dts+z[1], '-')
#plt.ylim([0.0036, 0.022])
#plt.savefig("HP_jamming_hplot2.png", dpi=400)
#plt.savefig("HP_jamming_hplot2.pdf")


plt.figure(2, figsize = (6,6))
for i in range(len(gs)):
  data = np.loadtxt(Efiles[i])
  mean = np.mean(data, axis=0, dtype=np.float32)
  y = (mean-gs[i])/Ls[i]
  erry = np.std(data, axis=0, dtype=np.float32)/Ls[i]
  print strs[i], ':'
  print erry/(mean-gs[i])*Ls[i]*100
  print '\n'
  plt.errorbar(dts, y, yerr=erry, fmt = stls[i], label = strs[i], linewidth=2, elinewidth=1.4, capsize=4.4)#, stls[i], markersize = 6, label = strs[i])
  #print 'sa:', data[-1,:]
  #print 'gs:', gs[i]/Ls[i]
  #plt.xlim([0,0.8])
  #plt.ylim(ylims[i])
plt.yscale('log')
plt.xscale('log')
plt.legend(loc = 'lower right')
plt.xlabel(" d$T$", fontsize = 16)
plt.ylabel("$E(T=0) - GS$", fontsize = 16)
plt.xlim([1.8e-6, 0.00105])
plt.ylim([0.0001,0.038])
plt.grid('on')
#plt.savefig('HP_jamming_scaling.pdf')
#plt.savefig('HP_jamming_scaling.png', dpi=400)

#data = np.loadtxt(files[1])
#plt.figure(2, figsize=(6,6))
plt.show()
