import numpy as np
from time import time
from tqdm import tqdm

# Solves problem, return False if it becomes unstable, True if it solve succesfully

# Solver function 
def BLSolver(ExpSolver, L, H, visc, dy, xyratio, uinf):

    dx = xyratio * dy

    x = np.arange(0,L,dx)
    y = np.arange(0,H,dy)

    y_size = np.shape(y)[0]
    x_size = np.shape(x)[0]

    # Velocity matrices
    v = np.zeros((x_size,y_size))
    u = np.zeros((x_size,y_size))

    # Boundary conditions
    ux0 = uinf
    uy0 = 0
    vx0 = 0
    vy0 = 0

    # Solver

    start_time = time()

    # Applying boundary conditions
    v[::,0] = np.ones(x_size) * vy0
    v[0] = np.ones(y_size) * vx0

    u[0] = np.ones(y_size) * ux0
    u[::,0] = np.ones(x_size) * uy0
    u[::,-1] = np.ones(x_size) * ux0

    if ExpSolver:
        
        # Construct matrices for problem vectorization
        A = (np.eye(y_size,k=-1) - 2*np.eye(y_size) + np.eye(y_size,k=1)) * visc / dy**2
        A = A[1:-1]

        B = ((np.eye(y_size,k=-1) - np.eye(y_size,k=1)))/(2*dy)
        B = B[1:-1]

        for i in tqdm(range(0,x_size-1), desc="Solving discrete field…", ascii=False, ncols=100):

            u[i+1,1:-1] = (dx /u[i,1:-1]) * (A @ u[i] - v[i,1:-1] * (B @ u[i])) + u[i,1:-1]

            for j in range(0,y_size-1):
                v[i+1,j+1] = xyratio * (u[i,j+1] - u[i+1,j+1]) + v[i+1,j]
            
            if max(u[i+1]) > 1.01:
                print('-------- INSTABILITY DETECTED, TERMINATING SOLVER --------')
                return False
            
        elapsed_time = time()-start_time

    elif ExpSolver == False:
        # Calculate some recurring values and construct some "primitive" matrices.    These are to achieve slightly faster solution times.
        a1 = visc / (2 * dy**2)
        a2 = 1 / (4 * dy)
        Mfv = np.eye(y_size,k=-1) * a2 - np.eye(y_size,k=1) * a2
        Mfu = np.eye(y_size) / dx
        Ck = (np.eye(y_size,k=-1) - 2*np.eye(y_size) + np.eye(y_size,k=1)) * a1

        for i in tqdm(range(0,x_size-1), desc="Solving discrete field…", ascii=False, ncols=100):

            # Prepare linear system
            K =  Mfv * (np.diag(v[i,0:-1],k=1) + np.diag(v[i,1::],k=-1)) + Ck

            B = Mfu * np.diag(u[i]) + K
            B[0,0]=1
            B[0,1]=0
            B[-1,-1]=1
            B[-1,-2]=0

            A = Mfu * np.diag(u[i]) - K
            A[0,0]=1
            A[0,1]=0
            A[-1,-1]=1
            A[-1,-2]=0

            # Solve
            u[i+1] = np.linalg.solve(A, B @ u[i])

            for j in range(0,y_size-1):
                v[i+1,j+1] = xyratio * (u[i,j+1] - u[i+1,j+1]) + v[i+1,j]
            
            if max(u[i+1]) > 1.01:
                print('-------- INSTABILITY DETECTED, TERMINATING SOLVER --------')
                return False

        elapsed_time = time()-start_time


    
    # Blasius Solution
    Rex = x * uinf /visc
    deltas = np.empty((3, x_size))
    deltas = np.array([[5], [1.72], [0.664]]) * np.concatenate(([0], (x[1::] / np.sqrt(Rex[1::]))))
    Cf = np.concatenate(([0], (0.664 / np.sqrt(Rex[1::]))))

    # Save data
    np.savez('BLData.npz', deltas = deltas, u = u, v = v, x = x, y = y, Cf = Cf, elapsed_time = elapsed_time)

    return True

