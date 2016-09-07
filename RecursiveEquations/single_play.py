# single plot

import numpy as np
import HNRG as hnrg
import matplotlib.pyplot as plt


hp = hnrg.hp6()

kl = hp.rg_all_kl(4, 0, 200)
#print kl

plt.plot(kl[:, 1], kl[:, 0], '-d')
plt.plot(kl[0, 1], kl[0, 0], 'rd', markersize=10)
plt.plot(kl[-1, 1], kl[-1, 0], 'cd', markersize=10)


plt.show()
