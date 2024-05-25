import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

pd_simi = pd.read_csv('./similarity_frag.csv')
pd_simi = pd_simi.drop(columns=['AP'])
#pd_simi = pd.DataFrame(pd_simi.values.T, index=pd_simi.columns, columns=pd_simi.index)
sns.set(color_codes=True)
#species = iris.pop("species")
#设置图片大小
g= sns.clustermap(pd_simi, fmt="d",cmap='YlGnBu',figsize=(12,10))
ax = g.ax_heatmap
label_y = ax.get_yticklabels()
plt.setp(label_y, rotation=360, horizontalalignment='left')

#设置图片名称，分辨率，并保存
plt.savefig('similarity.tiff',dpi = 600)
plt.show()