from __future__ import division
import os
import re
import numpy as np



def MakeInitialConditions(nparticles):
    vecs = np.zeros((nparticles,2,2))
    for vec in xrange(nparticles):
        r = ((12-5)*np.random.random())+5
        phi = 2*np.pi*np.random.random()
        phid = ((0.1-0.001)*np.random.random())+0.001
        vecs[vec] = [[r*np.cos(phi), r*np.sin(phi)],
                     [-r*np.sin(phi)*phid, r*np.cos(phi)*phid]]
    return(vecs)
        

def ElasticCollision(state,masses):
    raise ValueError("whoops!")


def G2xy(G,r,T):
    
    returnG = np.zeros_like(G)
    
    rd = G[0,0]
    rdd = G[0,1]
    Td = G[1,0]
    Tdd = G[1,1]
    
    xd = np.cos(T)*rd-r*np.sin(T)*Td
    xdd = np.cos(T)*(-r*(Td**2)+rdd)-np.sin(T)*(2*rd*Td+r*Tdd)
    yd = np.sin(T)*rd+np.cos(T)*r*Td
    ydd = 2*np.cos(T)*rd*Td+np.sin(T)*(-r*(Td**2)+rdd)+np.cos(T)*r*Tdd
    
    returnG[0] = [xd,yd]
    returnG[1] = [xdd,ydd]
    
    return(returnG)

def get_filenum(dir):
    filenums = [[0]]
    
    regex = re.compile(r'\d+')
    
    for filename in os.listdir(dir):
        
        if(filename[:5] == "nbody" and filename[-5:] == ".fits"):
            filenums.append([int(x) for x in regex.findall(filename[-7:])])
    return(1+int(max(filenums)[0]))



def xy2rad(state,forces):
    
    radstate = np.zeros([3,2])

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
        T = np.arctan(y/x)
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
        print state
        raise ValueError("NaN encountered")
    
    #velocities
    rd = (x*xd+y*yd)/r
    Td = (x*yd-y*xd)/r2
    
    #Accelerations
    rdd = (-4*(x*xd+y*yd)**2+4*r2*(xd**2+yd**2+x*xdd+y*ydd))/(4*(r2**(3/2)));
    Tdd = (y**2*(-xdd*y+2*xd*yd)-((x**2)*(xdd*y+2*xd*yd))+(x**3)*ydd+x*y*(2*(xd**2)-2*(yd**2)+y*ydd))/(r2**2)
    
    radstate[0] = [r,T]
    radstate[1] = [rd,Td]
    radstate[2] = [rdd-r*(Td**2),r*Tdd+2*rd*Td]
    
    return(radstate)

def a0(x):
        return(0)