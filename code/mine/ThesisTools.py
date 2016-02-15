from __future__ import division
import os
import re
import numpy as np

def SingleParticleNewtonianForce(state, masses, i, Nparticles, eps):
    forces = np.empty([Nparticles,2])
    
    x1 = state[i,0,0]
    y1 = state[i,0,1]
    m1 = masses[i]
    for particle_num in xrange(Nparticles):
        if particle_num != i:
            m2 = masses[particle_num]
            x2 = state[particle_num,0,0]
            y2 = state[particle_num,0,1]
            distance2 = ((x2-x1)**2+(y2-y1)**2)
            if distance2 < eps:
                jforce = m2/distance2
                jforcedir = [x2-x1,y2-y1]/np.sqrt(distance2)
                forces[particle_num] = jforce*jforcedir
    return(np.sum(forces,axis=0))

def G2xy(G,r,T):
    
    returnG = np.empty(G.shape)
    
    rd = G[0,0]
    rdd = G[0,1]
    Td = G[1,0]
    Tdd = G[1,1]
    
    xd = np.cos(T)*rd-r*np.sin(T)*Td
    xdd = np.cos(T)*(-r*Td**2+rdd)-np.sin(T)*(2*rd*Td+r*Tdd)
    yd = np.sin(T)*rd+np.cos(T)*r*Td
    ydd = 2*np.cos(T)*rd*Td+np.sin(T)*(-r*Td**2+rdd)+np.cos(T)*r*Tdd
    
    returnG[0] = [xd,yd]
    returnG[1] = [xdd,ydd]
    
    return(returnG)

def get_filenum(dir):
    filenums = [[0]]
    
    regex = re.compile(r'\d+.csv')
    
    for filename in os.listdir(dir):
        if filename[:5] == "nbody":
            filenums.append([int(x) for x in regex.findall(filename)])
    return(1+int(max(filenums)[0]))



def xy2rad(state,forces):
    
    radstate = np.empty([3,2])

    x = state[0,0]
    y = state[0,1]
    xd = state[1,0]
    yd = state[1,1]
    xdd = forces[0]
    ydd = forces[1]
    
    #positions
    r2 = x**2+y**2
    r = np.sqrt(r2)
    if x>0:
        T = np.arctan(y/x)+np.pi
    elif x<0:
        if y>=0:
            T = np.arctan(y/x)+np.pi
        elif y<0:
            T = np.arctan(y/x)-np.pi
    elif x==0:
        if y>0:
            T = np.pi/2
        elif y<0:
            T = -np.pi/2
        elif y==0:
            T = 0
    else:
        print x,y
        raise ValueError("NaN encountered")
    
    #velocities
    rd = (x*xd+y*yd)/r
    Td = (x*xd-y*yd)/r2
    
    #Accelerations
    rdd = (-4*(x*xd+y*yd)**2+4*r2*(xd**2+yd**2+x*xdd+y*ydd))/(4*(r2)**(3/2));
    Tdd = (y**2*(-xdd*y+2*xd*yd)-x**2*(xdd*y+2*xd*yd)+x**3*ydd+x*y*(2*xd**2-2*yd**2+y*ydd))/(r2**2)
    
    
    radstate[0] = [r,T]
    radstate[1] = [rd,Td]
    radstate[2] = [rdd-r*Td**2,r*Tdd+2*rd*Td]
    
    return(radstate)