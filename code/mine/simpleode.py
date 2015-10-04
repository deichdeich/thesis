import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.integrate import odeint

init_state = [.9,0]

def dstate_dt(state,t):
    x = state[0]
    xdot = state[1]
    dstate = np.zeros_like(state)
    
    dstate[0] = xdot
    #dstate[1] = -2 * x**3 * (1 - x * xdot - xdot) + (x**2 * xdot - x**3 * xdot)
    
    # simple precession test case
    # along with this x'', make init_state = [<1,~0]
    # dstate[1] = -x*(.5) + .5
    
    return dstate
    
t = np.arange(0.0,15*np.pi,0.001)

finalstate = odeint(dstate_dt,init_state,t)


r = 1/finalstate[:,0]
phi = t

x = r * np.cos(phi)
y = r * np.sin(phi)

fig = plt.figure()
plt.plot(x,y)
plt.show()