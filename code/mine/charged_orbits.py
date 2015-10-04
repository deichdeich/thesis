import numpy as np
from scipy.integrate import odeint

# [r,rdot,phi,phidot] 
init_state = [100.,1000.,np.pi,100.]


# define some constants
Jz = 1.
k = 1.
E = 1.
m = 1.
Q = 1.
L = 10.
    
def dstate_dt(state, t):    
    # unpack state vector
    r = state[0]
    rdot = state[1]
    phi = state[2]
    phidot = state[3]
    
    dstate = np.zeros_like(state)
    dstate[0] = rdot
    dstate[1] = ((r * phidot**2) + 1/(r**2) - (L * phidot)/(r**2))
    dstate[2] = phidot
    dstate[3] = L/(r**2) - rdot/r
    
    return dstate

t = np.arange(0.0,.2*np.pi,0.0001)

finalstate = odeint(dstate_dt, init_state, t)

# convert from spherical to cartesian
theta = np.pi/2
r = finalstate[:,0]
phi = finalstate[:,2]

x = r * np.sin(theta) * np.cos(phi)
y = r * np.sin(theta) * np.sin(phi)


import matplotlib.pyplot as plt
fig = plt.figure()
#plt.plot(x,y)
plt.scatter(x,y,c=np.log10(t),edgecolors='none',s=20)
plt.xlabel('X position')
plt.ylabel('Y position')
cb = plt.colorbar()
cb.set_label('log(time)')
plt.scatter(0,0,marker='x',s=100,color='black')

plt.show()
