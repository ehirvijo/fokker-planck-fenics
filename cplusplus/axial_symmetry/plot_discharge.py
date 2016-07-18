"""

Read in the discharge.dat file and plot the result of a run.

"""

import matplotlib.pyplot as plt
import numpy as np

data=np.loadtxt('discharge.dat',skiprows=0)

print np.shape(data)

j=data[:,3]
dt=data[1,1]-data[2,1]
dj=np.diff(j)/dt

f, axarr = plt.subplots(5,sharex=True)
axarr[0].plot(data[:,1],data[:,3],linewidth=2,label='J_{||}')
axarr[0].plot(data[:-1,1],dj,linewidth=2,label='J_{||}')
axarr[0].set_title("J_{||}")
axarr[1].plot(data[:,1],data[:,4],linewidth=2,label='E')
axarr[1].set_title("E")
axarr[2].plot(data[:,1],data[:,5],linewidth=2,label='Ti')
axarr[2].set_title("Ti")
axarr[3].plot(data[:,1],data[:,6],linewidth=2,label='ni')
axarr[3].set_title("ni")
axarr[4].plot(data[:,1],data[:,7],linewidth=2,label='ne')
axarr[4].set_title("ne")

plt.xlabel('t')
plt.show()


