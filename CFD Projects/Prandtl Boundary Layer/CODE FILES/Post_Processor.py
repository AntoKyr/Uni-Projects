import numpy as np
from matplotlib import pyplot as plt

with np.load('BLData.npz', allow_pickle=True) as BLData:
    deltas = BLData['deltas']
    u = BLData['u']
    v = BLData['v']
    x = BLData['x']
    y = BLData['y']
    Cfblas = BLData['Cf']

H = y[-1]
L = x[-1]

uinf = u[0,1]

x_size = np.shape(x)[0]
# Boundary layer thickness
ind1, full_range1 = np.where(u - 0.99*uinf > 0, u, np.inf).argmin(axis=1), np.arange(0,x_size)
ind2, full_range2 = np.where(u - 0.99*uinf < 0, u, -np.inf).argmax(axis=1), np.arange(0,x_size)

cfdelta0 = (0.99*uinf - u[full_range1,ind1]) * (y[ind2] - y[ind1]) / (u[full_range2,ind2] - u[full_range1,ind1]) + y[ind1]

# Displacement thickness 
cfdelta1 = np.zeros(x_size)
for i in range(0,x_size):
    cfdelta1[i] = np.trapz(1 - u[i] / uinf, y)

# Momentum thickness
cfdelta2 = np.zeros(x_size)
for i in range(0,x_size):
    cfdelta2[i] = np.trapz((u[i] / uinf) * (1 - u[i] / uinf), y)

misc=1.81*10**-5

Cfds = np.zeros(x_size)
for i in range(0,x_size):
    Cfds[i] = misc * (u[i,1]-u[i,0]) / (y[1]-y[0])

# Plots

# Boundary layer velocity distribution
plt.figure()
plt.imshow(np.rot90(u), cmap='gray', extent=(0,L,0,H), aspect='auto')

for i in range(0,np.shape(y)[0],3):
    plt.plot([x[-1], x[-1] + u[-1,i]], [y[i], y[i]], 'k')
plt.plot(u[-1] + x[-1], y, 'k')
plt.plot(x[-1] * np.ones(np.shape(y)[0]), y, 'k')

# Boundary layer thickness
plt.plot(x,deltas[0],'b',label = 'δ Blasius')
plt.plot(x,cfdelta0,'b--',label = 'δ cfd')

# Displacement thickness
plt.plot(x,deltas[1],'r',label = 'δ1 Blasius')
plt.plot(x,cfdelta1,'r--',label = 'δ1 cfd')

# Momentum thickness
plt.plot(x,deltas[2],'m',label = 'δ2 Blasius')
plt.plot(x,cfdelta2,'m--',label = 'δ2 cfd')

plt.axis([0,12,0,0.1])
plt.title('x-axis velocity distribution')
plt.xlabel('x [ m ]')
plt.ylabel('y [ m ]')
plt.legend()


plt.figure()
plt.plot(x,Cfds,'--',label='Cf cfd')
plt.plot(x,Cfblas, label='Cf Blasius')
plt.ylabel('Cf')
plt.xlabel('x [ m ]')
plt.legend()
plt.axis([0,11,0,0.01])
plt.grid()
plt.show()