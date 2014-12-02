import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

sns.set(style="darkgrid")

fig, axs = plt.subplots(3, figsize=(4, 12))
df = pd.read_csv(sys.argv[1])

axs[0].scatter(df.fisher, df.shuffled, s=4)
axs[0].set_xlim(0, 1)
axs[0].set_ylim(0, 1)
axs[0].set_xlabel('fisher p-value')
axs[0].set_ylabel('shuffled p-value')
axs[0].plot([0, 1], [0, 1], ls='--')

x = -np.log10(df.fisher)
y = -np.log10(df.shuffled)

m = int(max(x.max(), y.max())) + 1
axs[1].scatter(x, y, s=4)
axs[1].set_xlim(0, m)
axs[1].set_ylim(0, m)
axs[1].set_xlabel('-log10(fisher p-value)')
axs[1].set_ylabel('-log10(shuffled p-value)')
axs[1].plot([0, m], [0, m], ls='--')


x = -np.log10(1 - np.minimum(1-1e-6, df.fisher))
y = -np.log10(1 - np.minimum(1-1e-6, df.shuffled))

m = int(max(x.max(), y.max())) + 1
axs[2].scatter(x, y, s=4)
axs[2].set_xlim(0, m)
axs[2].set_ylim(0, m)
axs[2].set_xlabel('-log10(1 - fisher p-value)')
axs[2].set_ylabel('-log10(1 - shuffled p-value)')
axs[2].plot([0, m], [0, m], ls='--')

plt.tight_layout()
plt.savefig(sys.argv[1].replace('.txt', '') + '.png')

fig.show()
