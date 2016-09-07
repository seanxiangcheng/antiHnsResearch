
import numpy as np
import matplotlib.pyplot as plt

ks_1d = np.arange(10)[3:]
ks_pb = np.arange(16)[3:]
bad_1d = np.array([2,4,8,16,34,66,132])
bad_pb = np.array([2,4,8,16,33,66,132, 266, 532, 1064, 2128, 4257, 8514])

plt.figure(1, figsize = (6,6))
plt.semilogy(ks_pb,bad_pb,'-ko', linewidth=1.8, markersize=8, label='Periodic');
plt.semilogy(ks_1d, bad_1d, '-b>', linewidth=1.8, markersize=8, label='Open');
plt.xlim([np.amin(ks_pb)*0.95, np.amax(ks_pb)*1.01])
plt.ylim([np.amin(bad_pb)*0.95, np.amax(bad_pb)*1.05])
plt.grid('on')
plt.legend(loc='upper left')
plt.xlabel('k')
plt.ylabel('# of bad bonds')
plt.savefig('Bad_bonds_plot.pdf')
plt.show()

