import numpy as np
from matplotlib import pyplot as plt
from Solver import BLSolver


# Modes:
Mode =''
Mode = 'Stability'
#Mode = 'Convergence dy'
#Mode = 'Convergence dx/dy'
#Mode = 'Diagram'


if Mode == 'Stability':
    
    # Configure solver
    ExpSolver = False
    L = 5
    H = 0.1
    visc = 1.48 * (10**-5)
    uinf = 1

    Y_powers = -np.linspace(1.1,3,10)
    i=0

    # Empty arrays
    max_stable_ratios = np.zeros(np.shape(Y_powers)[0])
    

    for dy in 10**(Y_powers):

        xyedges = np.array([0.0,np.sqrt(2)],dtype=float)

        # Partition Method
        while xyedges[1]-xyedges[0] > 0.001:
            
            xycurr = np.average(xyedges)
            SolverBool = BLSolver(ExpSolver=ExpSolver, L=L, H=H, visc=visc, dy=dy, xyratio=xycurr, uinf=uinf)

            if SolverBool:
                xyedges[0] = xycurr
            else:
                xyedges[1] = xycurr

        max_stable_ratios[i] = xyedges[1]        
        i += 1

    # Diagram
    plt.plot(10**-Y_powers,max_stable_ratios)
    plt.title('Stability investigation')
    plt.ylabel('Maximum stable dx/dy')
    plt.xlabel('Number of nodes in 1m')
    plt.grid()
    plt.show()

#-----------------------------------------------------------------------------------------------------------------
elif Mode == 'Convergence dy':

    # Configure solver
    ExpSolver = False
    L = 11
    H = 0.15
    visc = 1.48 * (10**-5)
    uinf = 1
    xyratio = 0.7

    Y_powers = -np.linspace(1,3.2,10)
    i=0
    misc=1.81*10**-5

    # Empty arrays
    Cfds = np.zeros(np.shape(Y_powers)[0])
    
    for dy in 10**(Y_powers):

        BLSolver(ExpSolver=ExpSolver, L=L, H=H, visc=visc, dy=dy, xyratio=xyratio, uinf=uinf)

        with np.load('BLData.npz', allow_pickle=True) as BLData:
            u=BLData['u'][-1]

        # Cf compute
        Cfds[i] = misc * (u[1]-u[0]) / dy
        
        i += 1
    
    Rel = L * uinf /visc
    CfBlas = (0.664 / np.sqrt(Rel)) * np.ones(np.shape(Y_powers)[0])

    # Diagram
    plt.plot(10**-Y_powers, CfBlas, '--', label='Cf Blasius')
    plt.plot(10**-Y_powers, Cfds, label='Cf cfd')
    plt.title('Convergence investigation')
    plt.ylabel('Minimum Cf')
    plt.xlabel('Number of nodes in 1m')
    plt.grid()
    plt.show()

#-----------------------------------------------------------------------------------------------------------------
elif Mode == 'Convergence dx/dy':

    # Configure solver
    ExpSolver = True
    L = 11
    H = 0.15
    visc = 1.48 * (10**-5)
    uinf = 1
    dy = 0.005
    
    xyratio_powers = np.linspace(0,1.5,10)
    xyratios = 1 * 10**-(xyratio_powers)
    i=0
    misc=1.81*10**-5

    # Empty arrays
    Cfds = np.zeros(np.shape(xyratios)[0])
    
    for xyratio in xyratios:

        BLSolver(ExpSolver=ExpSolver, L=L, H=H, visc=visc, dy=dy, xyratio=xyratio, uinf=uinf)

        with np.load('BLData.npz', allow_pickle=True) as BLData:
            u = BLData['u'][-1]
    
        # Cf compute
        Cfds[i] = misc * (u[1]-u[0]) / dy
        
        i += 1
    
    Rel = L * uinf /visc
    CfBlas = (0.664 / np.sqrt(Rel)) * np.ones(np.shape(xyratios)[0])

    # Diagram
    plt.plot(1/xyratios, CfBlas, '--', label='Cf Blasius')
    plt.plot(1/xyratios, Cfds, label='Cf cfd')
    plt.title('Convergence investigation')
    plt.ylabel('Minimu Cf')
    plt.xlabel('dy / dx ratio')
    plt.grid()
    plt.show()

#-----------------------------------------------------------------------------------------------------------------
elif Mode == 'Diagram':

    # Configure solver
    ExpSolver = False
    L = 11
    H = 0.15
    visc = 1.48 * (10**-5)
    uinf = 1
    xyratio = 1
    dy = 0.001

    # Solve
    if BLSolver(ExpSolver=ExpSolver, L=L, H=H, visc=visc, dy=dy, xyratio=xyratio, uinf=uinf):

        # Run Post_Processor
        import Post_Processor
        Post_Processor