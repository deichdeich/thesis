import numpy as np
from scipy.integrate import odeint

# [r,theta,phi,rdot,thetadot,phidot] 
init_state = [100.,np.pi/2.,0.,1.,0.,1.]

# define some constants
Jz = 1.
k = 1.
E = 1.
m = 1.
Q = 1.
        
    
def dstate_dt(state, t):    
    # unpack state vector
    r = state[0]
    theta = state[1]
    phi = state[2]
    rdot = state[3]
    thetadot = state[4]
    phidot = state[5]
    
    dstate = np.zeros_like(state)
    dstate[0] = rdot
    dstate[1] = thetadot
    dstate[2] = phidot
    dstate[3] = (1./(2. * rdot)) * ((((2.*Jz**2.)/(m**2*r**3))*rdot) + k * (1./(r**2))*rdot + k * (((-1. * Jz**2)/(r**4)) * rdot + (k/(r**5))*rdot)) # r''
    dstate[4] = 0. # theta''
    dstate[5] = k/m * (1./(r**4) * rdot * phidot) - (2./r * rdot * phidot) # phi''

    return dstate

t = np.arange(0.0,50.,0.1)

finalstate = odeint(dstate_dt, init_state, t)

# convert from spherical to cartesian

r = finalstate[:,0]
theta = finalstate[:,1]
phi = finalstate[:,2]

x = r * sin(phi) * cos(theta)
y = r * sin(phi) * sin(theta)
z = r * cos(phi)

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot(x,y,z)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
plt.show()
