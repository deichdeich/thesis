from __future__ import division, print_function
import numpy128 as np
import time
import sys
import os
import ThesisTools as farts
from astropy.io import fits

class OrbitalSystem(object):
    def __init__(self,
                 Nparticles,
                 M,
                 a = farts.a0,
                 dt=0.01,
                 interaction = "ClassicalNBody",
                 masses = None,
                 init_state = None,
                 save_dir = "/Users/alexdeich/Dropbox/thesis/code/mine/nbody_output"):
        
        self.start_particles = Nparticles
        self.Nparticles = Nparticles
        self.interaction_type = interaction
        self.M = M
        self.dt = dt
        self.a = a
        self.event_horizon = 0
        self.cleanup = []
        if init_state is not None:
            self.use_state = np.copy(init_state)
        else:
            self.use_state = None
        self.save_dir = save_dir
        self.init_state = np.copy(init_state)
        self.state = np.copy(init_state)
        self.fname = ""
        
        if init_state is not None:
            if init_state.shape == (Nparticles,2,2):
                self.init_state = use_state
            else:
                raise ValueError("Initial state not the right shape: ",init_state.shape,
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
            rad = farts.xy2rad(kstate[particle],
                               farts.SingleParticleNewtonianForce(kstate,
                                                                  self.masses,
                                                                  particle,
                                                                  self.Nparticles,
                                                                  2))
                                                                  
                                                                  
        elif self.interaction_type == None:
            rad = farts.xy2rad(self.state[particle],(0,0))
            
        elif self.interaction_type == "ClassicalElastic":
            raise ValueError("ClassicalElastic collisions not implemented yet")
        
        
        r = rad[0,0]   
        phi = rad[0,1]
        
        
        
        f = np.array(([rad[0,0],rad[1,0]],[rad[0,1],rad[1,1]]))
        G=np.array([[f[0,1],
                     -(1/f[0,0]**4)*(self.a(t)**2-2*self.M*f[0,0]+(f[0,0]**2))*((self.M*(f[0,0]**4)*f[0,1]**2)/(self.a(t)**2-2*self.M*f[0,0]+(f[0,0]**2))**2-((f[0,0]**5)*(f[0,1]**2))/((self.a(t)**2)-2*self.M*f[0,0]+(f[0,0]**2))**2+self.M*(-1+self.a(t)*f[1,1])**2+(f[0,0]**3)*((f[0,1]**2)/((self.a(t)**2)-2*self.M*f[0,0]+(f[0,0]**2))-(f[1,1]**2)))+rad[2,0]],
                     [f[1,1],
                     -((2*f[0,1]*(self.a(t)*self.M+(-(self.a(t)**2)*self.M+f[0,0]**3)*f[1,1]))/(f[0,0]*(2*(self.a(t)**2)*self.M+(self.a(t)**2)*f[0,0]+(f[0,0]**3))))+rad[2,1]]])
        
        return(farts.G2xy(G,r,phi))
        
    def UpdateStateVectorRK4(self,t):
        self.event_horizon = self.get_event_horizon(t)
        new_state = np.ndarray((self.Nparticles,2,2))
        
        for particle in xrange(self.Nparticles):
            kstate = np.copy(self.state)
            #do RK4 shit
            k1 = self.dt * self.SingleParticleDerivativeVector(kstate,particle, t+self.dt)
            kstate[particle] = np.copy(self.state[particle]) + k1/2
            k2 = self.dt * self.SingleParticleDerivativeVector(kstate,particle,t+(self.dt/2))
            kstate[particle] = np.copy(self.state[particle]) + k2/2
            k3 = self.dt * self.SingleParticleDerivativeVector(kstate,particle,t+(self.dt/2))
            kstate[particle] = np.copy(self.state[particle]) + k3
            k4 = self.dt * self.SingleParticleDerivativeVector(kstate,particle,t+self.dt)
            new_state[particle] = np.copy(self.state[particle]) + (1/3)*(k1/2 + k2 + k3 + k4/2)
            
            r = np.sqrt(new_state[particle,0,0]**2+new_state[particle,0,1]**2)
            if(r > 999 or r < self.event_horizon):
                if particle not in self.cleanup:
                    self.cleanup.append(particle)
                    print("\n{} added to cleanup".format(particle))
        
        #Get rid of gobbled or ejected particles 
        if self.cleanup != []:
            for particle in self.cleanup:
                print("\n***particle {} shit the bed at step {}***".format(particle,int(t/self.dt)))
                self.remove_particle(particle)
                new_state = np.delete(new_state,particle,axis=0)
                print("***particle {} removed***".format(particle))
                self.cleanup.remove(particle)
                if self.cleanup == []:
                    print("***cleanup completed***")
        self.cleanup = []
        
        return(new_state)
    
    def MakeInitialConditions(self):
        self.cleanup = []
        vecs = np.zeros((self.start_particles,2,2))
        for vec in xrange(self.start_particles):
            r = ((12-5)*np.random.random())+5
            phi = 2*np.pi*np.random.random()
            phid = ((0.1-0.001)*np.random.random())+0.001
            vecs[vec] = [[r*np.cos(phi), r*np.sin(phi)],
                         [-r*np.sin(phi)*phid, r*np.cos(phi)*phid]]
        return(vecs)
        
        
    def TimeEvolve(self,nsteps,comments):
        self.cleanup = []
        t=0
        
        ##Get init_state
        if self.use_state == None:
            self.state = np.copy(self.MakeInitialConditions())
            self.init_state = np.copy(self.state)
            self.Nparticles = len(self.state)
        else:
            self.state = np.copy(self.use_state)
        print("Got initial conditions for {} particles".format(self.start_particles))
        
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
            avg = total_time/step
            perc = 100*((step+1)/nsteps)
            sys.stdout.write('\rFrame {} of {} completed ({}%).  Step: {}s, Total: {}s, Estimated time remaining: {}s. Nparticles: {}'.format(step+1,
                                                                       nsteps,
                                                                       '%0.1f'%perc,
                                                                       '%0.4f'%steptime,
                                                                       '%0.4f'%total_time,
                                                                       '%i'%(((avg * nsteps)+1)-total_time),
                                                                       '%i'%self.Nparticles))
            sys.stdout.flush()
        print("\nWriting to disk...")
        filenum = farts.get_filenum(self.save_dir)
        self.fname = "{}/nbody_{}_{}.fits".format(self.save_dir,self.start_particles,filenum)
        hdulist.writeto(self.fname,clobber=True)
        print("Data written at {}".format(self.fname))
    
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
    
    def trajplot(self,interval = "all",saveplot = False):
        framelist = fits.open(self.fname)
        
        if type(interval) == int:
            interval = [interval,interval+1]
        elif interval == "all":
            interval = [1,len(framelist)]
        x = []
        y = []
        rng = max(interval)-min(interval)
        if rng < 300:
            stepsize = 1
        elif rng < 20000:
            stepsize = 10
        else:
            stepsize = 500
        for frame in xrange(interval[0],interval[1],stepsize):
            x.append(float(framelist[frame].data[0][0]))
            y.append(float(framelist[frame].data[0][1]))
        
        import matplotlib.pyplot as plt
        plt.scatter(x,y, color="black",edgecolors= "None")
        plt.scatter(0,0,marker='x', color="black")
        circle1 = plt.Circle((0,0),radius = 2,color='r',fill=False)
        plt.xlim(-12,12)
        plt.ylim(-12,12)
        fig = plt.gcf()
        fig.gca().add_artist(circle1)
        plt.show()
    
    def movie(self):
        data = fits.open(self.fname)
        for i in xrange(1,len(data)):
            plotdata = np.rec.array([data[i].data["X"],data[i].data["Y"]],names=("x","y"))
            plt.figure()
            plt.scatter(plotdata["x"],plotdata["y"],alpha=0.3)
            plt.scatter(0,0,marker="x", color="black")
            circle1 = plt.Circle((0,0),radius = 2,color='r',fill=False)
            fig = plt.gcf()
            fig.gca().add_artist(circle1)
            plt.xlim(-30,30)
            plt.ylim(-30,30)
            fname = "/Users/alexdeich/30-particle-force-comparison/newton/frames/{}.png".format(i)
            plt.savefig(fname)
        
        os.system("ffmpeg -framerate 30 -i {}%d.png -c:v libx264 -r 30 -pix_fmt yuv420p out.mp4".format("/Users/alexdeich/30-particle-force-comparison/newton/frames/"))
         
    def get_event_horizon(self,t):
        return(self.M + np.sqrt(self.M**2 - self.a(t)**2))
    
    def densityplot(self,frame):
        raise ValueError("This is not implemented yet")
    
    def stats(self, interval, statistical_function):
        #If len(interval)==1: return histogram of that frame
        #Else: return plot of statistical_function over interval
        raise ValueError("This is not implemented yet")