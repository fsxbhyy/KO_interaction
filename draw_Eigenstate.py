from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import sys
from scipy import stats

np.set_printoptions(threshold=np.nan,suppress=True)

for i in range(len(sys.argv)):
    print sys.argv[i]

a=np.transpose(np.loadtxt(sys.argv[1],unpack=False))

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.scatter(a[0], a[1], (a[2]/a[0]), c='k',marker='o')

ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')

plt.show()
