import numpy as np
from matplotlib import pyplot as plt


def Solver(U, n, dt, ImplicitZeta, ImplicitPsi, Updates):
    # Problem Parameters
    L = 1
    v = 10**-5
    Re = abs(U*L/v)

    # Discretization
    dx = L/(n - 1)
    dy = L/(n - 1)

    # Relaxation
    omega = 1

    # Starting Conditions
    Zeta = np.zeros(n**2)
    Zeta[n-1::n] = -U/dy
    Zeta[n-2::n] = -U/(2*dy)
    Psi = np.zeros(n**2)

    # Constructing matrices for computation
    mat0 = np.eye(n**2)
    mat1 = np.eye(n**2,k=-n) + np.eye(n**2,k=n)
    mat2 = np.eye(n**2,k=-1) + np.eye(n**2,k=1)


    # Preparing Solver Operation
    NewZeta = Zeta

    ZetaError = np.array([Zeta[round(n**2/2)],10])
    ZetaCounter = 0
    ZetaConv = np.array([Zeta[round(n**2/2)]])
    PsiCounter = 0
    ConverVal = -1

    while (abs(Zeta[round(n**2/2)] - ZetaError[1]) > abs(Zeta[round(n**2/2)] - ZetaError[0]) * 10**-5) or (ZetaCounter<25):

        # Check for explosion or oscillation
        if (abs(Zeta[round(n**2/2)] - ZetaError[1])>15) or (ConverVal > abs(((Zeta[round(n**2/2)] - ZetaError[0]) * 10**-3) / (Zeta[round(n**2/2)] - ZetaError[1]))):
            print('-------------INSTABILITY DETECTED - TERMINATING SOLVER-------------')
            return False
        
        # Ready values for next iteration
        ConverVal = abs(((Zeta[round(n**2/2)] - ZetaError[0]) * 10**-3) / (Zeta[round(n**2/2)] - ZetaError[1]))
        ZetaError[1] = Zeta[round(n**2/2)]

        # Update on progress

        if Updates and (np.mod(ZetaCounter,20)==0):
            print(ConverVal)


        # Psi Calculation
        
        if ImplicitPsi:
            # linsys: A * Psi = B
            # Building A matrix and B vector
            Amat = mat2 / dy**2 + mat1 / dx**2 - 2 * mat0 * (1 / dx**2 + 1 / dy**2)
            Bvec = -Zeta

            # Clear all the boundary values
            BoundaryIndices = np.concatenate([np.arange(1,n-1), np.arange(n-1, n**2, n), np.arange(0, (n-1)*n+1, n), np.arange((n-1)*n+1, n**2-1)])

            Amat[BoundaryIndices,:] = 0
            Amat[:,BoundaryIndices] = 0
            for i in BoundaryIndices:
                Amat[i,i] = 1
            
            Bvec[BoundaryIndices] = 0

            # Solving linsys
            Psi = np.linalg.solve(a=Amat, b=Bvec)

        else:

            PsiError = np.array([Psi[round(n**2/2)],10])
            PsiCounter = 0

            while (abs(Psi[round(n**2/2)] - PsiError[1]) > abs(Psi[round(n**2/2)] - PsiError[0]) * 10**-4) or (PsiCounter<20):

                NewPsi = ((dx*dy)**2 * mat0@Zeta + dy**2 * mat1@Psi + dx**2 * mat2@Psi) / (2*(dx**2 + dy**2))
                PsiError[1] = Psi[round(n**2/2)]
                PsiCounter += 1

                # Psi = 0 at boundary
                NewPsi[1:n-1] = 0
                NewPsi[0 : (n-1)*n+1 : n] = 0
                NewPsi[n-1 : n**2 : n] = 0
                NewPsi[(n-1)*n+1 : n**2] = 0

                Psi = Psi + omega*(NewPsi - Psi)
    
        
        # Computing boundary Zeta values on walls
        NewZeta[1:n] = -2 * (Psi[n+1 : 2*n]-Psi[1:n]) / dx**2
        NewZeta[-n+1:] = 2 * (Psi[-n+1:] - Psi[-2*n+1:-n]) / dx**2
        NewZeta[0:-n+1:n] = -2 * (Psi[1:-n+2:n] - Psi[0:-n+1:n]) / dy**2

        # Computing boundary Zeta values on moving wall
        NewZeta[2*n-1:-n:n] = 2 * (-U/dy + (Psi[2*n-1:-n:n] - Psi[2*n-2:-n-1:n]) / dy**2)

        Zeta = NewZeta

        ZetaConv = np.concatenate((ZetaConv, np.array([Zeta[round(n**2/2)]])))
        ZetaCounter += 1

        # Predictor
        if ImplicitZeta:
            # linsys: A * Zeta = B
            # Building A matrix
            jmat = np.diag((-mat1@Psi / (8*dx*dy))[0:-1] - 1/(2 * Re * dy**2), k=1) + np.diag((mat1@Psi / (8*dx*dy))[1:] - 1/(2 * Re * dy**2), k=-1)
            imat = np.diag((mat2@Psi / (8*dx*dy))[0:-n] - 1/(2 * Re * dx**2), k=n) + np.diag((-mat2@Psi / (8*dx*dy))[n:] - 1/(2 * Re * dx**2), k=-n)
            centermat = mat0 * (1/dt + 1/(Re * dx**2) + 1/(Re * dy**2))
            Amat = imat+jmat+centermat

            # Building B vector
            Bvec = -(mat2@Psi * mat1@Zeta) / (8*dx*dy) + (mat1@Psi * mat2@Zeta) / (8*dx*dy) + ((mat1 - 2*mat0) @ Zeta / dx**2 + (mat2 - 2*mat0) @ Zeta / dy**2) / (2*Re) + Zeta / dt

            # Solving linsys
            NewZeta = np.linalg.solve(a=Amat, b=Bvec)

        else:
            NewZeta = dt * ((-mat2@Psi * mat1@Zeta + mat1@Zeta * mat2@Zeta) / (4*dy*dx) + ((mat1 - 2*mat0) @ Zeta / dx**2 + (mat2 - 2*mat0) @ Zeta / dy**2) / Re) + Zeta
    
    # Save data and kill
    np.savez('CFData.npz',n = n, Psi = Psi, Zeta = Zeta, U = U, L = L, ZetaCounter = ZetaCounter, ZetaConv = ZetaConv)
    return True