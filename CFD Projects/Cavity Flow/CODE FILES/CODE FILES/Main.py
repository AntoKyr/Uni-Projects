import numpy as np
from matplotlib import pyplot as plt
from CFDSolver import Solver


mode = 2

if mode == 0:
    # Density Convergence

    maxn = 35
    minn = 11
    PsiVal = np.zeros([round((maxn-minn)/8)+1])
    U = 500* 10**-5

    Counter = 0

    for n in range(minn,maxn+1,8):

        dt = abs(U * 10**3 * 3 / n)
        SolverBool = Solver(U = U, n = n, dt = dt, ImplicitZeta = True, ImplicitPsi = True, Updates = True)

        if SolverBool==False:
            break

        with np.load('CFData.npz', allow_pickle=True) as CFData:
            n = CFData['n']
            Psi = CFData['Psi']
        
        PsiVal[Counter] = Psi[round(n**2/2)] 
        Counter+=1
        print(Counter)
    
    # Diagram
    plt.figure()
    plt.plot(np.arange(minn,maxn+1,8),PsiVal[:], label = 'Re = 500')
    plt.title('Stream Function Convergence Re = 500')
    plt.ylabel('Center Stream Function Value')
    plt.xlabel('Number of nodes in 1m')
    plt.grid()

    plt.show()
            

elif mode == 1:
    # dt Stability 

    maxn = 31
    minn = 7
    nstep = 8
    max_stable_dt = np.zeros([round((maxn-minn)/nstep)+1,3])
    U = np.array([1,100,500])*10**-5

    for i in range(0,3):

        Counter = 0

        for n in range(minn,maxn+1,nstep):

            xyedges = np.array([0, 5],dtype=float)

            # Partition Method
            while xyedges[1]-xyedges[0] > 10**-4:
                
                xycurr = np.average(xyedges)
                SolverBool = Solver(U = U[i], n = n, dt = xycurr, ImplicitZeta = True, ImplicitPsi = True, Updates = True)
                print(10**-2/(xyedges[1]-xyedges[0]))

                if SolverBool:
                    xyedges[0] = xycurr
                else:
                    xyedges[1] = xycurr

            max_stable_dt[Counter,i] = xyedges[1]        
            Counter += 1
            print(Counter)

    # Diagram
    plt.plot(np.arange(minn,maxn+1,2),max_stable_dt[:,0], label = 'Re = 1')
    plt.plot(np.arange(minn,maxn+1,2),max_stable_dt[:,1], label = 'Re = 100')
    plt.plot(np.arange(minn,maxn+1,2),max_stable_dt[:,2], label = 'Re = 500')
    plt.title('Stability investigation')
    plt.ylabel('Maximum stable dt')
    plt.xlabel('Number of nodes in 1m')
    plt.legend()
    plt.grid()
    plt.show()

elif mode == 2:
    # 2D Field 
    n = 31
    U = 5000 * 10**-5
    dt = abs(U * 10**3 / (n * 15))
    Solver(U = U, n = n, dt = dt, ImplicitZeta = True, ImplicitPsi = True, Updates = True)
