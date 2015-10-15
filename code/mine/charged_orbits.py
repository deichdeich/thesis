import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# [r,rdot,phi,phidot] 
init_state = [2.,0.,0.,0.]


# define some constants
Jz = 1.
k = 1.
E = 1.
m = 1.
Q = 1.

        
def dstate_dt(state, t): 
    # unpack state vector
    
    if t<100:
        L = -0.005
    if t>100:
        L = 0.005
    
    #L = -0.005
    r = state[0]
    rdot = state[1]
    phi = state[2]
        
    dstate = np.zeros_like(state)
    dstate[0] = rdot
    dstate[1] = 1/r**3 - 1/r**2 + L * (-1/r**4 - L/(r**5))
    dstate[2] = 1/(r**2) + L/r**4
    return dstate
    
t = np.arange(0.1,200.,0.001)
    
finalstate = odeint(dstate_dt, init_state, t)
    
# convert from spherical to cartesian
theta = np.pi/2
r = finalstate[:,0]
phi = finalstate[:,2]
    
x = r * np.cos(phi)
y = r * np.sin(phi)
    
    
    
plt.figure()
plt.scatter(x,y,c=(t),s=10,edgecolors='none')
plt.xlabel('X position')
plt.ylabel('Y position')
cb = plt.colorbar()
cb.set_label('time')
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