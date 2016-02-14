from __future__ import division
import numpy as np
import ThesisTools as farts

class OrbitalSystem(object):
    def __init__(self,
                 Nparticles,
                 M,
                 a,
                 nsteps,
                 dt=0.001,
                 interaction_type = "ClassicalNBody",
                 masses = None,
                 use_state = None):
        
        self.Nparticles = Nparticles
        self.interaction_type = interaction_type
        self.M = M
        self.dt = dt
        self.nsteps = nsteps
        self.a = a
        self.event_horizon = 0
        self.cleanup = []
        
        if use_state == None:
            self.init_state = np.random.random((Nparticles,2,2))
        else:
            if use_state.shape == (Nparticles,2,2):
                self.init_state = use_state
            else:
                raise ValueError("Initial state not the right shape: ",use_state.shape,
                                 "\nShould be ",(Nparticles,2,2))
    
        if masses == None,
            self.masses = np.zeros(Nparticles)+0.0001
        else:
            if len(masses) == Nparticles:
                self.masses = masses
            else:
                raise ValueError("Mass list must be of length ", Nparticles,
                                 "\nGiven length:", len(masses))
    
    def SingleParticleDerivativeVector(self,kstate,particle,t):
        if self.interaction_type == "ClassicalNBody":
            rad = farts.xy2rad(state[particle],
                               farts.SingleParticleNewtonianForce(kstate,
                                                                  self.masses,
                                                                  self.Nparticles,
                                                                  2))
                                                                  
                                                                  
        elif self.interaction_type == None:
            rad = farts.xy2rad(state[particle],(0,0))
        elif self.interaction_type == "ClassicalElastic":
            raise ValueError("ClassicalElastic collisions not implemented yet")
        
        
        
        r = rad[0,0]
        
        if(r > 999 or r < self.event_horizon)):
            self.cleanup.append(particle)
           
        phi = rad[0,1]
        
        f = np.array((rad[0],rad[1]))
        
        G = np.array([[f[1],
                       -(1/f[0]**4)*(self.a(t)**2 - 2*self.M*f[0] + (f[0]**2))*((self.M*(f[0]**4)*f[1]**2)/(self.a(t)**2 - 2*self.M*f[0] + (f[0]**2))**2 - ((f[0]**5)*(f[1]**2))/((self.a(t)**2) - 2*self.M*f[0] + (f[0]**2))**2 + self.M*(-1 + self.a(t)*f[3])**2 +  (f[0]**3)*((f[1]**2)/((self.a(t)**2) - 2*self.M*f[0] + (f[0]**2)) - (f[3]**2))) + rself.ad[2, 0]],
                       [f[3],
                       -((2*f[1]*(self.a(t)*self.M + (-(self.a(t)**2)*self.M + f[0]**3)*f[3]))/(f[0]*(2*(self.a(t)**2)*self.M + (self.a(t)**2)*f[0] + (f[0]**3)))) + rself.ad[2, 1]]])
        
        
        return(G2xy(G,r,phi))
        
    def UpdateStateVectorRK4(self,t):
        self.event_horizon = self.M + np.sqrt(self.M**2 - self.a(t)**2)
        new_state = np.ndarray(self.state.shape)
        
        for particle in xrange(self.Nparticles):
            kstate = self.state
            #do RK4 shit
            k1 = dt * SingleParticleDerivativeVector(kstate,particle, t+self.dt)
            kstate[particle] = self.state[particle] + k1/2
            k2 = dt * SingleParticleDerivativeVector(kstate,particle,t+(self.dt/2))
            kstate[particle] = self.state[particle] + k2/2
            k3 = dt * SingleParticleDerivativeVector(kstate,particle,t+(self.dt/2))
            kstate[particle] = self.state[particle] + k3
            k4 = dt * SingleParticleDerivativeVector(kstate,particle,t+self.dt)
            new_state[particle] = self.state[particle] + (1/3)*(k1/2 + k2 + k3 + k4/2)
            
        #Get rid of gobbled or ejected particles 
        for particle in self.cleanup:
            self.remove_particle(particle) 
        self.cleanup = []
        
        return(new_state)
        
    
    def remove_particle(self,particle):
        
        self.state = np.delete(self.state,particle,axis=0)
        self.Nparticles -= 1
    
    def combine_particles(self,particle1,particle2)
        #get phase space of new particle
        
        #delete collided particles in state and mass
        #add new particle in state and mass
        
        #nparticles -= 1 
        
        mass1 = self.masses[particle1]
        mass2 = self.masses[particle2]
        
        xd1 = self.state[particle1,1,0]
        yd1 = self.state[particle1,1,1]
        
        xd2 = self.state[particle2,1,0]
        yd2 = self.state[particle2,1,1]
        
        totalmass = mass1+mass2
        
        newxd = ((mass1*xd1) + (m2*xd2))/totalmass
        newyd = ((mass1*yd1) + (m2*yd2))/totalmass
        
        self.state = np.delete(self.state,[particle1,particle2],axis=0)
        self.masses = np.delete(self.masses,[particle1,particle2],axis=0)
        
        new_phase_vector = [self.state[particle1,0,0],self.state[particle1,0,1],[newxd,newyd]]
        self.state = np.append(self.state,new_phase_vector)
        self.masses = np.append(self.masses,[mass1+mass2],axis=0)
        
        self.Nparticles -= 1