import numpy as np
import string

# Hard-coded potential function information
scale = 4.
zmin = -19.4442

def transitionPaths(numTrajPerSaddle, T=0.3, kB=0.001, gam=100, m=1):

    # Problem parameter initialization

    eVperJ = 6.24150934E18              # eV per J
    permol = 6.022E23                   # Number of particles per mole
    # Boltzmann constant in kJ / mol / K or eV / K
    kB = 1.381E-23                      # Boltzmann constant (J / K)
    en_units = "kj"
    kB = 1E-3*kB*permol if "j" in string.lower(en_units) else kB*eVperJ

    T = T                           # temperature
    gam = 50                        # inverse time constant
    m = 1                            # particle mass
    
    # Simulation parameters
    nsp = 2                          # # of energy landscape saddle pts.
    ntps = numTrajPerSaddle                          # total # trajectories to run per saddle point
    ntot = nsp*ntps                   # total # trajectories to run
    dt = 0.01                           # set timestep
    nt = 5000                        # max number of timesteps per trajectory
    tf = dt*nt                       # set final time
    brPt = np.zeros(2)                  # init vector of trajectory endpts
    cutoff = .1*kB*T                  # Set basin cutoff height above the minima
    print cutoff

    ixi = 1./(gam*m)                 # define as multiplicative factor
    D = kB*T*ixi                     # diffusion coefficient
    sq2D = np.sqrt(2*D*dt)                 # multiplicative factor for stochastic term
    
    # allocate space for trajectory vectors
    tempvec = np.zeros(2, dtype=(float,(nt,2)))
    traj = [] # squeeze(zeros([nsp,ntps], dtype=(float,(1,2))))
    
    # Grid discretization for plotting of potential function
    nx = 400                         # number of x spatial cells
    sx = 4.0                         # x-direction simulation size
    x = np.linspace(-sx, sx, num=2*nx, endpoint=True) # discretize grid with step size dx
    ny = 300                         # number of y spatial cells
    sy = 3.0                         # y-direction simulation size
    y = np.linspace(-sy, sy, num=2*ny, endpoint=True) # discretize grid with step size dx
        
    
    #---------------------------------------------------------------------------------------------------
    #
    # Set initial coordinates at random to one of the two saddle points. Initial values of y are 
    # determined by eye by observing the initial behavior of the low temperature trajectories.
    #
    #---------------------------------------------------------------------------------------------------
    tr0 = np.zeros([2, 2, nsp])
    y0 = np.array([0.3303, -0.4104])
    brPt = np.array([-1,-1])
    xdel = 0.03
    ydel = sq2D 
    # Generate starting point for each trajectory for each saddle point
    for i in range(0, nsp):
        tr0[0,:,i] = np.array([-xdel, y0[i]])
        tr0[1,:,i] = np.array([xdel, y0[i]])  
    
        #  + np.random.random()/3 - 1./6
    for i in range(0, nsp):
        for j in range(0, ntps):
            k = 0
            while k < 2:
                t = 0 
                tstep = 1
                yIso = ydel*np.random.randn(1)
                tr_init = np.array([tr0[k,0,i], yIso+tr0[k,1,i]])
                tr = tr_init
                while t < tf :
                    tempvec[k][tstep,:] = tr
                    tr = tr - dt*ixi*delH(tr) + sq2D*np.random.randn(2)
                    t += dt 
                    tstep = tstep + 1 
                    if inBasin(tr,cutoff):
                        print tr, H(tr)
                        if (tr[0] < 0 and k == 0) or (tr[0] > 0 and k == 1):
                            brPt[k] = tstep-1
                            break
                        else:
                            k -= 1
                            break
                k += 1
            reverseTemp = tempvec[0][1:brPt[0],:]
            reverseTemp = reverseTemp[::-1]
            temp = tempvec[1][1:brPt[1],:]
            traj.insert(i*nsp+j, np.vstack((reverseTemp, temp)))
    return traj


def inBasin( tr, threshold ):
    return ( H(tr) <= threshold )

def H(coord):
    x = coord[0]
    y = coord[1]
    return scale * ( -2*np.exp(-(5*x+2)**2/2-(5*y)**2/3) - 2*np.exp(-(5*x-2)**2/2-(5*y)**2/3) + \
    6*np.exp(-((5*x)**2+(5*y+0.2)**2)/2) - 6*np.exp(-((5*x)**2+(5*y)**2)/8) ) - zmin
        
def delH(coord):
    x = coord[0]
    y = coord[1]
    return np.array([dHdx(x,y), dHdy(x,y)])

# def dHdx(x,y):
#     return 5*scale * ( 2*(5*x+2)*np.exp(-(5*x+2)**2/2-(5*y)**2/3) + 2*(5*x-2)*np.exp(-(5*x-2)**2/2-(5*y)**2/3) - \
#     6*x*np.exp(-((5*x)**2+(5*y+0.2)**2)/2) + 1.5*x*np.exp(-((5*x)**2+(5*y)**2)/8) )

# def dHdy(x,y):
#     return 5*scale * ( (4/3)*y*np.exp(-(5*x+2)**2/2-(5*y)**2/3) + (4/3)*y*np.exp(-(5*x-2)**2/2-(5*y)**2/3) - \
#     6*(5*y+0.2)*np.exp(-((5*x)**2+(5*y+0.2)**2)/2) + 1.5*y*np.exp(-((5*x)**2+(5*y)**2)/8) )


def dHdx(x,y):
    return scale * ( 10*(5*x+2)*np.exp(-(5*x+2)**2/2-(5*y)**2/3) + 10*(5*x-2)*np.exp(-(5*x-2)**2/2-(5*y)**2/3) - \
    150*x*np.exp(-((5*x)**2+(5*y+0.2)**2)/2) + 37.5*x*np.exp(-((5*x)**2+(5*y)**2)/8) )

def dHdy(x,y):
    return scale * ( (100/3)*y*np.exp(-(5*x+2)**2/2-(5*y)**2/3) + (100/3)*y*np.exp(-(5*x-2)**2/2-(5*y)**2/3) - \
    30*(5*y+0.2)*np.exp(-((5*x)**2+(5*y+0.2)**2)/2) + 37.5*y*np.exp(-((5*x)**2+(5*y)**2)/8) )
