import numpy as np
from matplotlib import pyplot as plt
import csv 
from tqdm import tqdm


### Solver
# GridDef is the number of nodes in the grid
# L is the length of the nozzle
# Afunction is the A = f(x) function
# FlType is the Flow type,
# 0 for subsonic flow
# 1 for supersonic isentropic
# 2 for supersonic with sock
# Eq type are the types of equations. True for Conservative False for Non-Conservative
def CFDSolver(A, GridDef, Steps, Courant, EqType, FlType, pe, Cx): 

    
    x = np.linspace(0,1,GridDef)
    rho = np.ones((Steps,GridDef))
    rho[0] = rho[0] - 0.3146 * x 
    T = np.ones((Steps,GridDef)) 
    T[0] = pe / rho[0]
    V = np.ones((Steps,GridDef)) * 0.1
    V[0] = (V[0] + 1.09 * x) * np.sqrt(T[0])
    drho = np.ones((Steps-1,GridDef)) 
    dT = np.ones((Steps-1,GridDef)) 
    dV = np.ones((Steps-1,GridDef)) 
    Alog = np.log(A)

    dx = 1 / (GridDef - 1)
    g = 1.4
    R = 8.314

    if FlType == 0:
        Me2 = 2*(pe ** -((g - 1) / g) - 1) / (g - 1)
        As = A[-1] / ( ( 1 / Me2) * ((2 / (g + 1)) * (1 + 0.5 * (g - 1) * Me2)) ** ((g + 1) / (g - 1)))
        A = A / As

    dt = np.empty((Steps-1))

    # Computation Matrix

    B = np.eye(GridDef, k=1) - np.eye(GridDef, k=-1)
    B[0,1] = 0

    B = (1 / (2 * dx) ) * B

    Sn = np.eye(GridDef, k=1) - 2 * np.eye(GridDef, k=0) + np.eye(GridDef, k=-1)
    Sn[0,0] = 0; Sn[0,1] = 0
    Sd = np.eye(GridDef, k=1) + 2 * np.eye(GridDef, k=0) + np.eye(GridDef, k=-1)

    # Boundary Matrices for calculating edge values
    Vb = np.eye((GridDef))
    Vb[0,0] = 0; Vb[0,1] = 1; Vb[0,2] = 1; Vb[0,3] = -1
    Vb[-1,-1] = 0; Vb[-1,-2] = 1; Vb[-1,-3] = 1; Vb[-1,-4] = -1
    rhob = np.eye((GridDef))
    rhob[-1,-1] = 0; rhob[-1,-2] = 1; rhob[-1,-3] = 1; rhob[-1,-4] = -1
    Tb = np.eye((GridDef))
    Tb[-1,-1] = 0; Tb[-1,-2] = 1; Tb[-1,-3] = 1; Tb[-1,-4] = -1

    # Equation types
    if EqType == True:
        # Con Solution
        for i in tqdm(np.arange(1,Steps), desc="Solving mess…", ascii=False, ncols=100):
            # Computing Variables
            U11 = rho[i-1] * A 
            U21 = rho[i-1] * A * V[i-1]
            U31 = rho[i-1] * (T[i-1] / (g - 1) + 0.5 * g * (V[i-1] ** 2)) * A

            F11 = U21
            F21 = (U21 ** 2) / U11 + ((g - 1) / g) * (U31 - (g / 2) * ((U21 ** 2) / U11))
            F31 = g * (U21 * U31 / U11) - 0.5 * g * (g-1) * (U21 ** 3 / U11 ** 2)
            J1 = ((g - 1) / g) * (U31 - (g / 2) * U21 ** 2 / U11) * (B @ Alog)

            # Predictor
            dU11 = - B @ F11
            dU21 = - B @ F21 + J1
            dU31 = - B @ F31

            # Compute dt
            dt[i-1] = Courant * np.nanmin(dx / (V[i-1] + np.sqrt(T[i-1])))

            # Estimator
            S = Cx * np.abs(Sn @ (rho[i-1] * T[i-1])) / (Sd @ (rho[i-1] * T[i-1]))
            U12 = dU11 * dt[i-1] + U11 + S * (Sn @ U11)
            U22 = dU21 * dt[i-1] + U21 + S * (Sn @ U21)
            U32 = dU31 * dt[i-1] + U31 + S * (Sn @ U31)

            F12 = U22
            F22 = (U22 ** 2) / U12 + ((g - 1) / g) * (U32 - (g / 2) * ((U22 ** 2) / U11))
            F32 = g * (U22 * U32 / U12) - 0.5 * g * (g-1) * (U22 ** 3 / U12 ** 2)
            J2 = ((g - 1) / g) * (U32 - (g / 2) * U22 ** 2 / U12) * (B @ Alog)

            dU12 = - B @ F12
            dU22 = - B @ F22 + J2
            dU32 = - B @ F32

            # Corrector
            dU1 = (dU11 + dU12) / 2
            dU2 = (dU21 + dU22) / 2
            dU3 = (dU31 + dU32) / 2

            U12 = dU1 * dt[i-1] + U11 + S * (Sn @ U12)
            U22 = dU2 * dt[i-1] + U21 + S * (Sn @ U22) 
            U32 = dU3 * dt[i-1] + U31 + S * (Sn @ U32) 

            rho[i] = rhob @ (U12 / A)
            T[i] = Tb @ ((g - 1) * (U32 / U12 - 0.5 * g * V[i] ** 2))
            V[i] = Vb @ (U22 / U12)

            # Subsonic Correction
            if FlType != 1:
                rho[i,-1] = np.sqrt(rho[i,-1] * pe / T[i,-1])
                T[i,-1] = pe / rho[i,-1]        

            drho[i-1] = (rho[i] - rho[i-1]) / dt[i-1]
            dV[i-1] = (V[i] - V[i-1]) / dt[i-1]
            dT[i-1] = (T[i] - T[i-1]) / dt[i-1]

    elif EqType == False:
        # Non-Con Solution
        for i in tqdm(np.arange(1,Steps), desc="Solving mess…", ascii=False, ncols=100):
            # Predictor
            drho1 = -rho[i-1] * (B @ V[i-1]) - rho[i-1] * V[i-1] * (B @ Alog) - V[i-1] * (B @ rho[i-1])
            dV1 = -V[i-1] * (B @ V[i-1]) - (1 / g) * (B @ T[i-1] + (T[i-1] / rho[i-1]) * (B @ rho[i-1]))
            dT1 = -V[i-1] * (B @ T[i-1]) - (g - 1) * T[i-1] * ((B @ V[i-1]) + V[i-1] * (B @ Alog))

            # Compute dt
            dt[i-1] = Courant * np.nanmin(dx / (V[i-1] + np.sqrt(T[i-1]))) 

            # Estimator
            S = Cx * np.abs(Sn @ (rho[i-1] * T[i-1])) / (Sd @ (rho[i-1] * T[i-1]))
            rho1 = rhob @ (rho[i-1] + drho1 * dt[i-1] + S * (Sn @ rho[i-1]))
            V1 = Vb @ (V[i-1] + dV1 * dt[i-1] + S * (Sn @ V[i-1]))
            T1 = Tb @ (T[i-1] + dT1 * dt[i-1] + S * (Sn @ T[i-1])) 

            drho2 = - rho1 * (B @ V1) - rho1 * V1 * (B @ Alog) - V1 * (B @ rho1)
            dV2 = - V1 * (B @ V1) - ( 1 / g) * (B @ T1 + (T1 / rho1) * (B @ rho1))
            dT2 = - V1 * (B @ T1) - (g - 1) * T1 * ((B @ V1) + V1 * (B @ Alog))

            # Corrector
            drho[i-1] = (drho1 + drho2) / 2
            dV[i-1] = (dV1 + dV2) / 2
            dT[i-1] = (dT1 + dT2) / 2

            S = Cx * np.abs(Sn @ (rho1 * T1)) / (Sd @ (rho1 * T1))
            rho[i] = rhob @ (rho[i-1] + drho[i-1] * dt[i-1] + S * (Sn @ rho1))
            T[i] = Tb @ (T[i-1] + dT[i-1] * dt[i-1] + S * (Sn @ T1))
            V[i] = Vb @ (V[i-1] + dV[i-1] * dt[i-1] + S * (Sn @ V1))

            # Subsonic Correction
            if FlType != 1:
                rho[i,-1] = np.sqrt(rho[i,-1] * pe / T[i,-1])
                T[i,-1] = pe / rho[i,-1]



    OutputData = [rho, V, T, drho, dV, dT]
    return OutputData


# Analytical solution
def AnalSolver(A, GridDef, pe):
    pass




# Pre and Post processing
L = 3
GridDef = 160
def Afunction(duct):
    x = np.linspace(0,L,GridDef)
    if duct == 0:
        A = 1 + 2.2 * (x - 1.5) ** 2
    else:
        A1 = 1 + 2.2 * (x[x<1.5] - 1.5) ** 2
        A2 = 1 + 0.2223 * (x[x>=1.5] - 1.5) ** 2
        A = np.concatenate((A1,A2))
    return A
A = Afunction(0)
Courant = 0.4
EqType = True
FlType = 1
Steps = 50
Cx = 0
pe =0.93

OutputData = CFDSolver(A=A, GridDef=GridDef, Courant=Courant, EqType=EqType, FlType=FlType, Steps=Steps, pe=pe, Cx=Cx)


rho = OutputData[0]
#print(rho[-10:-1])
plt.figure()
plt.plot(np.arange(0,Steps), rho[::,15])
plt.grid()
plt.title('rho')
plt.xlabel('Steps')
plt.ylabel('V / a0')
plt.figure()
plt.grid()
plt.plot(np.linspace(0,1,GridDef), rho[-1])
plt.show()







