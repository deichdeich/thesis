from __future__ import division, print_function
import numpy as np
import time
import sys
import ThesisTools as farts
from astropy.io import fits

class OrbitalSystem(object):
    def __init__(self,
                 Nparticles,
                 M,
                 a = None,
                 dt=0.001,
                 interaction_type = "ClassicalNBody",
                 masses = None,
                 use_state = None,
                 save_dir = "/Users/alexdeich/Dropbox/thesis/code/mine/nbody_output"):
        
        self.start_particles = Nparticles
        self.Nparticles = Nparticles
        self.interaction_type = interaction_type
        self.M = M
        self.dt = dt
        self.a = a
        self.event_horizon = 0
        self.cleanup = []
        self.use_state = use_state
        self.save_dir = save_dir
        
        if use_state is not None:
            if use_state.shape == (Nparticles,2,2):
                self.init_state = use_state
            else:
                raise ValueError("Initial state not the right shape: ",use_state.shape,
                                 "\nShould be ",(Nparticles,2,2))
    
        if masses == None:
            self.masses = np.zeros(Nparticles)+0.0001
        else:
            if len(masses) == Nparticles:
                self.masses = masses
            else:
                raise ValueError("Mass list must be of length ", Nparticles,
                                 "\nGiven length:", len(masses))
    
    def SingleParticleDerivativeVector(self,kstate,particle,t):
        if self.interaction_type == "ClassicalNBody":
            rad = farts.xy2rad(self.state[particle],
                               farts.SingleParticleNewtonianForce(kstate,
                                                                  self.masses,
                                                                  particle,
                                                                  self.Nparticles,
                                                                  2))
                                                                  
                                                                  
        elif self.interaction_type == None:
            rad = farts.xy2rad(state[particle],(0,0))
        elif self.interaction_type == "ClassicalElastic":
            raise ValueError("ClassicalElastic collisions not implemented yet")
        
        
        
        r = rad[0,0]
        
        if(r > 999 or r < self.event_horizon):
            self.cleanup.append(particle)
           
        phi = rad[0,1]
        
        f = np.array((rad[0],rad[1]))
        G = np.array([[f[0,1],
                       -(1/f[0,0]**4)*(self.a(t)**2 - 2*self.M*f[0,0] + (f[0,0]**2))*((self.M*(f[0,0]**4)*f[0,1]**2)/(self.a(t)**2 - 2*self.M*f[0,0] + (f[0,0]**2))**2 - ((f[0,0]**5)*(f[0,1]**2))/((self.a(t)**2) - 2*self.M*f[0,0] + (f[0,0]**2))**2 + self.M*(-1 + self.a(t)*f[1,1])**2 +  (f[0,0]**3)*((f[0,1]**2)/((self.a(t)**2) - 2*self.M*f[0,0] + (f[0,0]**2)) - (f[1,1]**2))) + rad[2, 0]],
                       [f[1,1],
                       -((2*f[0,1]*(self.a(t)*self.M + (-(self.a(t)**2)*self.M + f[0,0]**3)*f[1,1]))/(f[0,0]*(2*(self.a(t)**2)*self.M + (self.a(t)**2)*f[0,0] + (f[0,0]**3)))) + rad[2, 1]]])
        
        
        return(farts.G2xy(G,r,phi))
        
    def UpdateStateVectorRK4(self,t):
        self.event_horizon = self.M + np.sqrt(self.M**2 - self.a(t)**2)
        new_state = np.ndarray(self.state.shape)
        
        for particle in xrange(self.Nparticles):
            kstate = self.state
            #do RK4 shit
            k1 = self.dt * self.SingleParticleDerivativeVector(kstate,particle, t+self.dt)
            kstate[particle] = self.state[particle] + k1/2
            k2 = self.dt * self.SingleParticleDerivativeVector(kstate,particle,t+(self.dt/2))
            kstate[particle] = self.state[particle] + k2/2
            k3 = self.dt * self.SingleParticleDerivativeVector(kstate,particle,t+(self.dt/2))
            kstate[particle] = self.state[particle] + k3
            k4 = self.dt * self.SingleParticleDerivativeVector(kstate,particle,t+self.dt)
            new_state[particle] = self.state[particle] + (1/3)*(k1/2 + k2 + k3 + k4/2)
            
        #Get rid of gobbled or ejected particles 
        for particle in self.cleanup:
            self.remove_particle(particle) 
        self.cleanup = []
        
        return(new_state)
    
    def MakeInitialConditions(self):
        vecs = np.ndarray((self.Nparticles,2,2))
        for vec in xrange(self.Nparticles):
            r = ((10-5)*np.random.random())+5
            phi = 2*np.pi*np.random.random()
            phid = ((1.5-0.5)*np.random.random())+0.5
            vecs[vec] = [[r*np.cos(phi), r*np.sin(phi)],
                         [-r*np.sin(phi)*phid, r*np.cos(phi)*phid]]
        return(vecs)
        
        
    def TimeEvolve(self,nsteps,comments):
        t=0
        ##Get init_state
        if self.use_state == None:
            self.state = self.MakeInitialConditions()
        
        print("Generated initial conditions")
        
        
        primary = self.get_header(nsteps,comments)
        frame0 = self.get_hdu()
        
        hdulist = fits.HDUList([primary,frame0])
        total_time = 0
        for step in xrange(1, nsteps):
            stepstart = time.time()
            self.state = self.UpdateStateVectorRK4(t)
            framen = self.get_hdu()
            hdulist.append(framen)
            t += self.dt
            end = time.time()
            steptime = end-stepstart
            total_time += steptime
            sys.stdout.write("\rFrame {} of {} completed ({}%).  Last frame took {} seconds, {} seconds total".format(step+1,
                                                                                                    nsteps,
                                                                                                    ((step+1)/nsteps)*100,
                                                                                                    steptime,
                                                                                                    total_time))
        
        filenum = farts.get_filenum(self.save_dir)
        fname = "{}/nbody_{}_{}.fits".format(self.save_dir,self.start_particles,filenum)
        hdulist.writeto(fname,clobber=True)
    
    def get_header(self,nsteps,comments=""):
        
        prihdr = fits.Header()
        prihdr["NPARTS"] = self.Nparticles
        prihdr["INTERACT"] = self.interaction_type
        prihdr["DT"] = self.dt
        prihdr["BH_M"] = self.M
        if self.a:
            prihdr["SPINFUNC"] = True
        else:
            prihdr["SPINFUNC"] = False
        prihdr["NSTEPS"] = nsteps
        prihdr["COMMENTS"] = comments
        prihdu = fits.PrimaryHDU(header=prihdr)
    
        return(prihdu)
        
    def get_hdu(self):
        state_transpose = self.state.T
        frame = fits.BinTableHDU.from_columns([fits.Column(name='X',format='20A',array = state_transpose[0][0]),
                                                fits.Column(name='Y',format='20A',array = state_transpose[1][0]),
                                                fits.Column(name='Xd',format='20A',array = state_transpose[0][1]),
                                                fits.Column(name='Yd',format='20A',array = state_transpose[1][1])])
        return(frame)
        
    def remove_particle(self,particle):
        
        self.state = np.delete(self.state,particle,axis=0)
        self.Nparticles -= 1
    
    def combine_particles(self,particle1,particle2):
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
        self.state = np.append(self.state,new_phase_vector,axis=0)
        self.masses = np.append(self.masses,totalmass,axis=0)
        
        self.Nparticles -= 1
        
    def trajplot(self,interval):
        raise ValueError("This is not implemented yet")
    
    def densityplot(self,frame):
        raise ValueError("This is not implemented yet")
    
    def stats(self, interval, statistical_function):
        #If len(interval)==1: return histogram of that frame
        #Else: return plot of statistical_function over interval
        raise ValueError("This is not implemented yet")