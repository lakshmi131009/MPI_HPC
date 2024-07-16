import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


df = pd.read_csv('~/output.csv',header=None)
df1=df.head(90)
df2=df[91:181]
df3=df[181:271]
df4=df.tail(90)

df1["cat"] = np.tile(np.repeat(["send", "pack", "derived"], 5), 6)
df2["cat"] = np.tile(np.repeat(["send", "pack", "derived"], 5), 6)
df3["cat"] = np.tile(np.repeat(["send", "pack", "derived"], 5), 6)
df4["cat"] = np.tile(np.repeat(["send", "pack", "derived"], 5), 6)



df1["N"] = np.repeat([16,32,64,128,256,512], 15)
df2["N"] = np.repeat([16,32,64,128,256,512], 15)
df3["N"] = np.repeat([16,32,64,128,256,512], 15)
df4["N"] = np.repeat([16,32,64,128,256,512], 15)

ax = sns.boxplot(data=df1, x="N", y=0, hue="cat", zorder=0.9, dodge=False)
sns.pointplot(data=df1, x="N", y=0, hue="cat", estimator=np.median, ci=None, markers="None", ax=ax)

ax.set_title('p=16')
ax.set_xlabel('N^2')
ax.set_ylabel('time in seconds')

plt.savefig('plot16.jpg')


ax = sns.boxplot(data=df2, x="N", y=0, hue="cat", zorder=0.9, dodge=False)
sns.pointplot(data=df2, x="N", y=0, hue="cat", estimator=np.median, ci=None, markers="None", ax=ax)

ax.set_title('p=36')
ax.set_xlabel('N^2')
ax.set_ylabel('time in seconds')

plt.savefig('plot36.jpg')

ax = sns.boxplot(data=df3, x="N", y=0, hue="cat", zorder=0.9, dodge=False)
sns.pointplot(data=df3, x="N", y=0, hue="cat", estimator=np.median, ci=None, markers="None", ax=ax)
ax.set_title('p=49')
ax.set_xlabel('N^2')
ax.set_ylabel('time in seconds')

plt.savefig('plot49.jpg')


ax = sns.boxplot(data=df4, x="N", y=0, hue="cat", zorder=0.9, dodge=False)
sns.pointplot(data=df4, x="N", y=0, hue="cat", estimator=np.median, ci=None, markers="None", ax=ax)
ax.set_title('p=64')
ax.set_xlabel('N^2')
ax.set_ylabel('time in seconds')
plt.savefig('plot64.jpg')

