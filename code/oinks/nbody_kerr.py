"""
Orbit Integrator for Numerical Kerr Solutions (OINKS)
"It's a real resource hog!"

Author: Alex Deich
Date: February, 2016
"""

from __future__ import division, print_function
import numpy128 as np
import time
import matplotlib.pyplot as plt
import sys
import os
import ThesisTools as farts
import scipy.stats
from astropy.io import fits

class OrbitalSystem(object):
    def __init__(self,
                 Nparticles,
                 M=1,
                 a = None,
                 dt=0.01,
                 interaction = False,
                 collisions = 'elastic',
                 masses = None,
                 init_state = None,
                 save_dir = "./output"):
        
        self.start_particles = Nparticles
        self.Nparticles = Nparticles
        self.interaction = interaction
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
        self.collision_dict = {}
        self.collisions=collisions
        self.collided_particles = np.array([])
        self.fname_list = []
        self.dirname = ""
        self.plotdata = 0
        self.skip_mkdir=False
        self.nsteps = 0
        

        if self.Nparticles == None and self.init_state is not None:
            self.Nparticles = len(self.init_state)
            
        if self.a == None:
            self.a = farts.a0
        
        if init_state is not None:
            if init_state.shape == (Nparticles,2,2):
                self.init_state = init_state
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
        
        if self.interaction == False:
            rad = farts.xy2rad(self.state[particle],(0,0))
        else:
            rad = farts.xy2rad(kstate[particle],
                               self.SingleParticleNewtonianForce(particle,
                                                                 5,
                                                                 0.001))
        
        
        r = rad[0,0]   
        phi = rad[0,1]
        
        if(r > 999 or r < self.event_horizon):
            if particle not in self.cleanup:
                self.cleanup.append(particle)        
        
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
            #sys.stdout.write("\r{}".format(particle))
            #sys.stdout.flush()
            kstate = np.copy(self.state)
            #do RK4 shit
            #is copying the arrays necessary? this might be causing massive slowdowns
            k1 = self.dt * self.SingleParticleDerivativeVector(kstate,particle, t)
            kstate[particle] = np.copy(self.state[particle]) + k1/2
            k2 = self.dt * self.SingleParticleDerivativeVector(kstate,particle,t+(self.dt/2))
            kstate[particle] = np.copy(self.state[particle]) + k2/2
            k3 = self.dt * self.SingleParticleDerivativeVector(kstate,particle,t+(self.dt/2))
            kstate[particle] = np.copy(self.state[particle]) + k3
            k4 = self.dt * self.SingleParticleDerivativeVector(kstate,particle,t+self.dt)
            new_state[particle] = np.copy(self.state[particle]) + (1/3)*(k1/2 + k2 + k3 + k4/2)
            
        
        #Get rid of gobbled or ejected particles 
        
        if self.cleanup != []:
            new_state = np.delete(new_state,self.cleanup,axis=0)
            self.remove_particle(particle)
            for particle in self.cleanup:
                print("\n***particle {} shit the bed at step {}***".format(particle,int(t/self.dt)))
                print("***particle {} removed***".format(particle))
                self.cleanup.remove(particle)
                if self.cleanup == []:
                    print("***cleanup completed***")
        self.cleanup = []
        
        
        big_particle_list = []
        big_vector_list = []
        new_masses_list = []
        
        

        #Cleanup collided particles
        if self.collisions == "inelastic":
            if self.collision_dict != {}:       
                for key in self.collision_dict:
                    particles = self.collision_dict[key][1]
                    m = 0
                    x = self.collision_dict[key][0][0]
                    y = self.collision_dict[key][0][1]
                    mvx = 0
                    mvy = 0
                    for particle in particles:
                        m += self.masses[particle]
                    for particle in particles:
                        mvx += (self.state[particle][1][0]*self.masses[particle])
                        mvy += (self.state[particle][1][1]*self.masses[particle])
                        big_particle_list.append(particle)
                    new_masses_list.append(m)
                new_particle = [[x,y],[mvx/m,mvy/m]]
                big_vector_list.append(new_particle)

                self.masses = np.delete(self.masses,big_particle_list,axis=0)
                self.masses = np.append(self.masses,np.array(new_masses_list),axis = 0)
                new_state = np.delete(new_state,big_particle_list,axis=0)
                new_state = np.append(new_state,np.array(big_vector_list),axis = 0)
                self.Nparticles = len(new_state)
                self.collision_dict = {}
                self.collided_particles = np.array([])
                
                
                """
                For each particle in a given collision,
                get the distance to every other particle
                get the closest particle
                evaluate that particle in an elastic collision
                up the state vector accordingly
                """
                
            elif self.collisions == "elastic":
                if self.collision_dict != {}:
                    checked_particles = []
                    for key in self.collision_dict:
                        particles = self.collision_dict[key][1]
                        for i in particles:
                            if i not in checked_particles:
                                checked_particles.append(i)
                                distances = {}
                                for j in particles:
                                    if j != i and j not in checked_particles:
                                        checked_particles.append(j)
                                        distances[j] = np.sqrt((self.state[i][0][0]-self.state[j][0][0])**2+(self.state[i][0][1]+self.state[j][0][1])**2) # adds an entry in the distances dictionary with the particle number as the key and the value as the distance to the particle being evaluated
                                    k = min(distances, key = distances.get) # get the particle number of the closest particle
                                    particle1 = self.state[i]
                                    m1 = masses[particle1]
                                    x1 = particle1[0][0]
                                    x1d = particle1[1][1]
                                    y1 = particle1[0][1]
                                    y1d = particle1[1][1]
                                    particle2 = self.state[k]
                                    m2 = masses[particle1]
                                    x2 = particle2[0][0]
                                    x2d = particle2[1][1]
                                    y2 = particle2[0][1]
                                    y2d = particle2[1][1]
                                    newp1d = [(x1d*(m1-m2) + 2*m2*x2d)/(m1+m2),(y1d*(m1-m2) + 2*m2*y2d)/(m1+m2)]
                                    newp2d = [(x2d*(m1-m2) + 2*m2*x1d)/(m1+m2),(y2d*(m1-m2) + 2*m2*y1d)/(m1+m2)]
                                    newp1 = [[x1,y1],newp1d]
                                    newp2 = [[x2,y2],newp2d]
                                    new_state[i] = newp1
                                    new_state[k] = newp2
                self.collision_dict = {}
        
        new_state = new_state[new_state[:,0,0].argsort()] # arsort() returns the indices that sort the array.  Here, I'm sorting by x position     
        return(new_state)
    
    
    def TimeEvolve(self,nsteps,comments,write=True):
        self.cleanup = []
        self.nsteps = nsteps
        t=0
        
        ##Get init_state
        if self.use_state == None:
            self.cleanup = []
            self.state = np.copy(farts.make_initial_conditions2(self.start_particles,self.M,self.a(0)))
            self.init_state = np.copy(self.state)
            self.Nparticles = len(self.state)
        else:
            self.state = np.copy(self.init_state)
        print("Got initial conditions for {} particles".format(self.start_particles))
        
        primary = self.get_header(nsteps,comments)
        frame0 = self.get_hdu()
        
        filenum = farts.get_filenum(self.save_dir,self.start_particles)
        self.dirname = "{}/nbody_{}_{}".format(self.save_dir,self.start_particles,filenum)
        os.mkdir(self.dirname)
        os.mkdir("{}/data".format(self.dirname))
        hdulist = fits.HDUList([primary,frame0])
        total_time = 0
        savenums = nsteps/1000
            
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
            if step%1000 == 0:
            
                if write == True:
                    print("\nWriting to disk...")
                    fname = "{}/data/{}.fits".format(self.dirname,step)
                    hdulist.writeto(fname,clobber=True)
                    print("Frames {} - {} written at {}".format(step-1000,step,fname))
                    hdulist = fits.HDUList([primary])
                    self.fname_list.append(fname)
            
        if len(hdulist)!=1:
            print("\nWriting to disk...")
            fname = "{}/data/{}.fits".format(self.dirname,step)
            hdulist.writeto(fname,clobber=True)
            print("Frames {} written at {}".format(step+1,fname))
            self.fname_list.append(fname)
                        
    def get_header(self,nsteps,comments=""):
        prihdr = fits.Header()
        prihdr["NPARTS"] = self.Nparticles
        prihdr["INTERACT"] = self.interaction
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
        frame = fits.BinTableHDU.from_columns([fits.Column(name='X',format='E',array = state_transpose[0][0]),
                                                fits.Column(name='Y',format='E',array = state_transpose[1][0]),
                                                fits.Column(name='Xd',format='E',array = state_transpose[0][1]),
                                                fits.Column(name='Yd',format='E',array = state_transpose[1][1]),
                                                fits.Column(name='MASS',format='E',array = self.masses)])
        return(frame)
        
    def remove_particle(self,particle):
        self.state = np.delete(self.state,self.cleanup,axis=0)
        for key in self.collision_dict:
            if particle in self.collision_dict[key][1]:
                self.collision_dict[key][1].remove(particle)
        self.Nparticles = len(self.state)
    
    
    def SingleParticleNewtonianForce(self, i, soi_radius, collision_radius):
        
        forces = np.zeros([self.Nparticles,2])
        
        x1 = self.state[i,0,0]
        y1 = self.state[i,0,1]
        
        m1 = self.masses[i]
        collided1 = 0
        
        if self.interaction=='ClassicalNBody':
            xmin = x1-soi_radius
            xmax = x1+soi_radius
            ymin = y1-soi_radius
            ymax = y1+soi_radius
            sliced_indices_x = np.where(np.logical_and(self.state[:,0,0]>xmin,self.state[:,0,0]<xmax))
            sliced_indices_y = np.where(np.logical_and(self.state[:,0,1]>ymin,self.state[:,0,1]<ymax))
            sliced_indices = np.intersect1d(sliced_indices_x,sliced_indices_y)
            sliced_arr = self.state[sliced_indices]
        
            no_i_indices = np.where(np.logical_and(sliced_arr[:,0,0] == x1,sliced_arr[:,0,1]==y1))
            sliced_arr = np.delete(sliced_arr,no_i_indices,axis=0)
            sliced_masses = self.masses[sliced_indices]
            sliced_masses = np.delete(sliced_masses,no_i_indices,axis=0)
        
            distances2 = np.array(np.transpose(np.matrix((sliced_arr[:,0,0]-x1)**2+(sliced_arr[:,0,1]-y1)**2)))
            jforce = np.array(np.transpose(np.matrix(sliced_masses)))/distances2
            jforcedir = np.array(np.transpose(np.matrix([sliced_arr[:,0,0]-x1,sliced_arr[:,0,1]-y1])))/np.sqrt(distances2)
            forces = jforce*jforcedir
        
            """
            Collision calculations
        
            if i has not been noted as having collided:
            
            
        
            if i not in self.collided_particles:
                #print(i)
                self.collision_dict[i] = [(x1,y1),0]
                x_collision_min = x1-collision_radius
                x_collision_max = x1+collision_radius
                y_collision_min = y1-collision_radius
                y_collision_max = y1+collision_radius
        
                collision_indices_x = np.where(np.logical_and(self.state[:,0,0]>x_collision_min,
                                                              self.state[:,0,0]<x_collision_max))
                #print(collision_indices_x)
                collision_indices_y = np.where(np.logical_and(self.state[:,0,1]>y_collision_min,
                                                              self.state[:,0,1]<y_collision_max))
                #print(collision_indices_y)
                collision_indices = np.intersect1d(collision_indices_x, collision_indices_y)
                #print(collision_indices)
                if i in self.collided_particles:
                    collision_indices = np.delete(collision_indices,np.where(collision_indices==i))
                #print(collision_indices)
                collision_arr = self.state[collision_indices]
        
                new_particles = np.setdiff1d(collision_indices,self.collided_particles)
                self.collided_particles = np.append(self.collided_particles,new_particles)
                self.collision_dict[i][1] = new_particles
                #if len(new_particles)>0:
                #   print(self.collision_dict)
                #  print(self.collided_particles)

        
            """
            """
            Any way to get rid of these for loops and/or replace them with numpy?
            """
            for particle_num in sliced_indices:
                if particle_num != i:
                    collided2 = 0
                    m2 = self.masses[particle_num]
                    x2 = self.state[particle_num,0,0]
                    y2 = self.state[particle_num,0,1]
                    distance2 = ((x2-x1)**2+(y2-y1)**2)
                    #if self.interaction == "ClassicalNBody":
                            #jforce = m2/distance2
                            #jforcedir = [x2-x1,y2-y1]/np.sqrt(distance2)
                            #forces[particle_num] = jforce*jforcedir
                    """
                    if the collision dict is empty, put in the two particles being evaluated.
                
                    if it's not empty, first search to see if the main particle is in there
                    if the main particle is in there, move on to the secondary
                    if the secondary particle is not in there, put it in the entry with the main particle
                    """
                
                    if self.collisions != False and i not in self.cleanup and particle_num not in self.cleanup:
                        if distance2 < collision_radius:
                            if self.collision_dict == {}:
                                self.collision_dict[0] = [(x1,y1),[i,particle_num]]
                            else:
                                i_val = self.dict_check(self.collision_dict, i)
                                keynum = max(self.collision_dict)+1
                                if i_val == -1:
                                    self.collision_dict[keynum] = [(x1,y1),[i]]
                                    i_val = keynum
                                elif self.dict_check(self.collision_dict, particle_num) == -1:
                                    self.collision_dict[i_val][1].append(particle_num)

            return(np.sum(forces,axis=0))
        
        elif self.interaction == False:
            xmin = x1-collision_radius
            xmax = x1+collision_radius
            ymin = y1-collision_radius
            ymax = y1+collision_radius
            sliced_indices_x = np.where(np.logical_and(self.state[:,0,0]>xmin,self.state[:,0,0]<xmax))
            sliced_indices_y = np.where(np.logical_and(self.state[:,0,1]>ymin,self.state[:,0,1]<ymax))
            sliced_indices = np.intersect1d(sliced_indices_x,sliced_indices_y)
            sliced_arr = self.state[sliced_indices]
        
            no_i_indices = np.where(np.logical_and(sliced_arr[:,0,0] == x1,sliced_arr[:,0,1]==y1))
            sliced_indices = np.delete(sliced_indices, no_i_indices)
            
            for particle_num in sliced_indices:
                if self.collisions != False and i not in self.cleanup and particle_num not in self.cleanup:
                        if self.collision_dict == {}:
                            self.collision_dict[0] = [(x1,y1),[i,particle_num]]
                        else:
                            i_val = self.dict_check(self.collision_dict, i)
                            keynum = max(self.collision_dict)+1
                            if i_val == -1:
                                self.collision_dict[keynum] = [(x1,y1),[i]]
                                i_val = keynum
                            elif self.dict_check(self.collision_dict, particle_num) == -1:
                                self.collision_dict[i_val][1].append(particle_num)
            return((0,0))
    
    
    
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
        elif rng < 21000:
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
    
    def dict_check(self, some_dict, thing):
        what_key= -1
        for key in some_dict:
            if thing in some_dict[key][1]:
                what_key = key
        return what_key
    
    def get_event_horizon(self,t):
        return(self.M + np.sqrt(self.M**2 - self.a(t)**2))
        
        
        
        
	#####  Plotting Tools  #####
    
    def movie(self,type=None,function = "count"):
        if self.skip_mkdir == False:
            os.mkdir("{}/movies".format(self.dirname))
            os.mkdir("{}/movies/spatial".format(self.dirname))
            os.mkdir("{}/movies/spatial/frames".format(self.dirname))
            os.mkdir("{}/movies/heatmap".format(self.dirname))
            os.mkdir("{}/movies/heatmap/frames".format(self.dirname))
        if type == 'spatial':
            self.make_spatial()
        elif type == 'heatmap':
            self.make_heatmap(function)
        elif type == None:
            self.make_heatmap(function)
            self.make_spatial()
    
    def make_heatmap(self,function):
        print("\n")
        print("Creating animated heatmap...")
        plt.figure()
        n=1
        for fname in self.fname_list:
            wholedata = fits.open(fname)
            for i in xrange(1,len(wholedata)):
            	data = wholedata[i].data
                bins = scipy.stats.binned_statistic_2d(data['X'],
                                                       data['Y'],
                                                       data['MASS'],
                                                       statistic=function,
                                                       bins=50,
                                                       range=[[-30,30],[-30,30]])
				
                x,y = np.where(~np.isnan(bins[0]))
                plt.scatter(x,y,c=bins[0][np.where(~np.isnan(bins[0]))],
                                                            marker='s',
                                                            edgecolors='none',
                                                            s=200)
                t = n*self.dt
                plt.xlim(0,48)
                plt.ylim(0,48)
                circle1 = plt.Circle((24,24), radius = self.get_event_horizon(t), color = 'r', fill = False)
                fig = plt.gcf()
                fig.gca().add_artist(circle1)
                fig.gca().axes.get_xaxis().set_visible(False)
                fig.gca().axes.get_yaxis().set_visible(False)
                plt.axes().set_aspect('equal')
                plt.title('Nparticles: {}'.format(len(data['X'])))
                cb = plt.colorbar()
                plt.clim(0,50)
                cb.set_label('Number density')
                plt.scatter(24,24,marker="x", color="white")
                fname = "{}/movies/heatmap/frames/{}.png".format(self.dirname,n)
                plt.savefig(fname)
                plt.clf()
                sys.stdout.write('\rFrame {} completed'.format(n))
                sys.stdout.flush()
                n+=1  
            wholedata.close()
        os.system("ffmpeg -framerate 300 -i {}/movies/heatmap/frames/%d.png -c:v libx264 -r 30 -pix_fmt yuv420p {}/movies/heatmap/out.mp4".format(self.dirname,self.dirname))

    def make_spatial(self):
        print("\n")
        print("Creating position space movie...")
        plt.figure()
        n=1
        for fname in self.fname_list:
            wholedata = fits.open(fname)
            for i in xrange(1,len(wholedata)):
                data = wholedata[i].data
                plt.scatter(data['X'],data['Y'],c=data['MASS'])
                cb = plt.colorbar()
                plt.clim(0,50)
                cb.set_label('Particle mass')
                plt.scatter(0,0,marker="x", color="black")
                t = n*self.dt
                circle1 = plt.Circle((0,0), radius = self.get_event_horizon(t), color = 'r', fill = False)
                fig = plt.gcf()
                fig.gca().add_artist(circle1)
                fig.gca().axes.get_xaxis().set_visible(False)
                fig.gca().axes.get_yaxis().set_visible(False)
                plt.xlim(-30,30)
                plt.ylim(-30,30)
                plt.axes().set_aspect('equal')
                plt.title('Nparticles: {}'.format(len(data['X'])))
                fname = "{}/movies/spatial/frames/{}.png".format(self.dirname,n)
                plt.savefig(fname)
                plt.clf()
                sys.stdout.write('\rFrame {} completed'.format(n))
                sys.stdout.flush()
                n+=1
            wholedata.close()
        print('\n')
        os.system("ffmpeg -framerate 300 -i {}/movies/spatial/frames/%d.png -c:v libx264 -r 30 -pix_fmt yuv420p {}/movies/spatial/out.mp4".format(self.dirname,self.dirname))

    def plot_n(self):
        ns = np.ndarray(len(self.nsteps))
        i=0
        for fname in self.fname_list:
            data = fits.open(fname)
            ns[i] = len(data[1].data)-1
            i+=1
        ns = ns/self.init_particles
        plt.plot(ns)
        plt.set_xlabel('Timestep')
        plt.set_ylabel('\r$n_step/N$')
    
    def plot_rd(self):

    	r_arr = np.ndarray(self.nsteps)
    	rd_arr = np.ndarray(self.nsteps)
    	n=1
    	for fname in self.fname_list:
    		wholedata = fits.open(fname)
    		for i in xrange(1,len(wholedata)):
    			data = wholedata[i].data
    			r = np.mean(np.sqrt(data['X']**2+data['Y']**2))
    			r_arr[n] = r
    			rd = (1/r)*(data['X']*data['Xd']+data['Y']*data['Yd'])
    			rd_arr[n] = np.mean(rd)
    			n+=1
    		print(fname)
    	r_arr = r_arr/r_arr[0]
    	plt.plot(r_arr,linewidth=2,color='black')
    	plt.xlabel(r'Timestep')
    	plt.ylabel(r'$\frac{<r>}{<r_0>}$',fontsize=16)
    	plt.show()
    	plt.plot(rd_arr,linewidth=2,color='black')
    	plt.xlabel(r'Timestep')
    	plt.ylabel(r'$<\dot{r}>$',fontsize=16)
    	plt.show()