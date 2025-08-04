import numpy as np
from matplotlib import pyplot as plt


with np.load('CFData.npz', allow_pickle=True) as CFData:
    n = CFData['n']
    Psi = CFData['Psi']
    Zeta = CFData['Zeta']
    U = CFData['U']
    L = CFData['L']
    ZetaCounter = CFData['ZetaCounter']
    ZetaConv = CFData['ZetaConv']


dx = L/(n-1)
dy = L/(n-1)
# Velocity Computation

Vx = (np.eye(n**2,k=1) - np.eye(n**2,k=-1))@Psi / (2*dy)
Vy = -(np.eye(n**2,k=-n) - np.eye(n**2,k=n))@Psi / (2*dx)

Vx[1:n-1] = 0
Vx[0:-n+1:n] = 0
Vx[2*n-1:-n: n] = U
Vx[-n:] = 0
Vy[1:n-1] = 0
Vy[0:-n+1:n] = 0
Vy[2*n-1:-n: n] = 0
Vy[-n:] = 0

# Post-Processor

PsiField = np.rot90(np.reshape(Psi,(n,n)),k=-1)
ZetaField = np.rot90(np.reshape(Zeta,(n,n)),k=-1)
VxField = np.rot90(np.reshape(Vx,(n,n)),k=-1)
VyField = np.rot90(np.reshape(Vy,(n,n)),k=-1)
VField = np.sqrt(VyField**2 + VxField**2)

SizeFact = 1/(n*U)

y = np.linspace(0,L,n)
x = np.linspace(0,L,n)
X,Y = np.meshgrid(x,y)

# Zeta Convergence
plt.figure()
plt.plot(np.arange(0,ZetaCounter),ZetaConv[1:])
plt.grid()
plt.title('Vorticity Value Convergence')
plt.xlabel('Number of Iterations')
plt.ylabel('Central Zeta Value')

# Stream Contours
plt.figure()
plt.contour(X, Y, PsiField,20,cmap='jet')
plt.grid()
plt.title('Stream Function Contours')
plt.xlabel('x [m]')
plt.ylabel('y [m]')
plt.colorbar()

# Vorticity Contours
plt.figure()
plt.contour(X, Y, ZetaField,20,cmap='jet')
plt.grid()
plt.title('Vorticity Function Contours')
plt.xlabel('x [m]')
plt.ylabel('y [m]')
plt.colorbar()

# Velocity Field
plt.figure()
plt.streamplot(X, Y, VxField, VyField, color=VField, cmap='jet')
plt.colorbar(label='[m/s]')
plt.grid()
plt.title('Velocity Field')
plt.xlabel('x [m]')
plt.ylabel('y [m]')

plt.figure()
plt.quiver(X-VxField*SizeFact/2, Y-VyField*SizeFact/2, VxField*SizeFact, VyField*SizeFact)
plt.grid()
plt.title('Velocity Vector Field')
plt.xlabel('x [m]')
plt.ylabel('y [m]')

plt.show()