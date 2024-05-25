import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

pd_simi = pd.read_csv('./similarity_frag.csv')
pd_simi = pd_simi.drop(columns=['AP'])

print(pd_simi.shape[1],type(pd_simi.shape[1]))
print(pd_simi.shape[0],type(pd_simi.shape[0]))

plt.figure(figsize=(3, 2), dpi=600) 
np_simi = np.array(pd_simi)
np_simi = np_simi.flatten()
sns.distplot(np_simi,norm_hist=False)
plt.tick_params('both',width=2,labelsize=8)
plt.savefig('./similarity_distribution.tiff')
plt.show()
