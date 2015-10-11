import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# [r,rdot,phi,phidot] 
init_state = [.02,88.,0.001,-90.]


# define some constants
Jz = 1.
k = 1.
E = 1.
m = 1.
Q = 1.

        
def dstate_dt(state, t): 
    # unpack state vector
    L = 1.
    r = state[0]+0.000001
    rdot = state[1]
    phi = state[2]
    phidot = state[3]
        
    dstate = np.zeros_like(state)
    dstate[0] = rdot
    dstate[1] = ((r * phidot**2) - 1/(r**2) + (L * phidot)/(r**2))
    dstate[2] = phidot
    dstate[3] = -(L/(r**2) - rdot/r)
        
    return dstate
    
t = np.arange(0.1,1.,0.00001)
    
finalstate = odeint(dstate_dt, init_state, t)
    
# convert from spherical to cartesian
theta = np.pi/2
r = finalstate[:,0]
phi = finalstate[:,2]
    
x = r * np.cos(phi)
y = r * np.sin(phi)
    
    
    
plt.figure()
plt.scatter(x,y,c=np.log10(t),s=10,edgecolors='none')
plt.xlabel('X position')
plt.ylabel('Y position')
cb = plt.colorbar()
cb.set_label('log(time)')
plt.scatter(0,0,marker='x',s=100,color='black')
plt.show()
"""
plt.title('init_state:{}'.format(init_state))
plt.scatter(x,y,c=np.log10(t),s=10,edgecolors='none')
plt.xlabel('X position')
plt.ylabel('Y position')
cb = plt.colorbar()
cb.set_label('log(time)')
plt.scatter(0,0,marker='x',s=100,color='black')
plt.savefig('{}_spatial.png'.format(init_state))
    
plt.figure()
plt.title('init_state:{}'.format(init_state))
plt.plot(t,r)
plt.xlim(9.9,10)
plt.xlabel('time')
plt.ylabel('distance from center')
plt.savefig('{}_radius.png'.format(init_state))
"""