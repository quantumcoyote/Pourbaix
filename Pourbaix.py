import subprocess
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate

#
# Change the range to change pH. Remember 0 to 15 is ph 0 to p
#

for i in range(0,15,1):
    pH=float(i)
    subprocess.call("python ReactionMatrix.py "+str(10**(-1*float(i))), shell=True)
    P=1.0
    for T in range(270,370,10):
        subprocess.call("python Thermochemistry.py "+str(T)+" "+str(P)+" "+str(pH), shell=True)
        subprocess.call("python StandardState.py "+str(T)+" "+str(P)+" "+str(pH), shell=True)

subprocess.call("mv momentinertia.dat scratch/.", shell=True)

x=[]
y=[]
z=[]
plot=open("scratch/plot.dat","r")
for line in plot:
    tmp=line.split()
    y.append(float(tmp[0]))
    x.append(float(tmp[1]))
    z.append(float(tmp[2]))

x=np.asarray(x)
y=np.asarray(y)
z=np.asarray(z)

xi, yi = np.linspace(x.min(), x.max(), 100), np.linspace(y.min(), y.max(), 100)
xi, yi = np.meshgrid(xi, yi)
rbf = scipy.interpolate.Rbf(x, y, z, function='linear')
zi = rbf(xi, yi)

plt.imshow(zi, vmin=z.min(), vmax=-z.min(),origin='lower',extent=[x.min(), x.max(), y.min(), y.max()],aspect='auto',cmap='coolwarm')
#plt.imshow(zi, vmin=z.min(), vmax=z.max(), origin='lower', extent=[x.min(), x.max(), y.min(), y.max()],aspect='auto', cmap='coolwarm')

plt.xlabel('pH')
plt.ylabel('Temperature (K)')
plt.colorbar()
plt.savefig("plot.png")
subprocess.call("rm scratch/*", shell=True)


