from __future__ import division
import os
import numpy as np
from math import sqrt

def elastic_collision(state,masses):
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

def get_filenum(dir,npart):
    nums = [0]
    for filename in os.listdir(dir):
        if(filename[:5] == "nbody"):
            pieces = filename.split("_")
            if int(pieces[-2]) == npart:
                nums.append(int(pieces[-1]))
    return(max(nums)+1)


def make_initial_conditions(nparticles):
    
    # this is a stable-ish circular orbit: [[[7,0],[0,.35]]]
    vecs = np.zeros((nparticles,2,2))
    for vec in xrange(nparticles):
        r = np.random.normal(10,.5)
        phi = (2*np.pi)*np.random.random()
        phid = abs(np.random.normal(0.05,0.15))
        vecs[vec] = [[r*np.cos(phi), r*np.sin(phi)],
                     [-r*np.sin(phi)*phid, r*np.cos(phi)*phid]]
    return(vecs)


def make_initial_conditions2(nparticles,bh_m,a0):
    vecs = np.zeros((nparticles,2,2))
    for vec in xrange(nparticles):
        r_init = np.nan
        while np.isnan(r_init):
        	energy = np.random.normal(5,.05)
        	jz = np.random.normal(7,1)
        	r_init = (1/(8 * bh_m)) * (4*(a0**2) - (a0**2) * energy**2 + jz**2 + np.sqrt((-(a0**2) * (-4 + energy**2) + jz**2)**2 - 48 * (-(a0**2) * energy + jz)**2 * bh_m**2))
        r_init += np.random.normal(0.5,0.25)
        phi = 2*np.pi*np.random.random(1)
        rd = 0
        phid = (pow(a0,5)*bh_m + 2*pow(a0,3)*bh_m*r_init*(-2*bh_m + r_init) + a0*bh_m*pow(r_init,2)*pow(-2*bh_m + r_init,2) - np.sqrt(-(pow(r_init,3)*pow(pow(a0,2) + r_init*(-2*bh_m + r_init),2)*(pow(a0,4)*bh_m*(-1 + pow(rd,2)) - pow(a0,2)*r_init*(2*bh_m*r_init + pow(r_init,2)*pow(rd,2) + pow(bh_m,2)*(-4 + pow(rd,2))) + bh_m*pow(r_init,2)*(-4*pow(bh_m,2) + 4*bh_m*r_init + pow(r_init,2)*(-1 + pow(rd,2)))))))/((pow(a0,2)*bh_m - pow(r_init,3))*pow(pow(a0,2) + r_init*(-2*bh_m + r_init),2))
        vecs[vec] = [[r_init*np.cos(phi), r_init*np.sin(phi)],
                    [-r_init * np.sin(phi) * phid, np.cos(phi)* r_init * phid]]
    return(vecs)

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