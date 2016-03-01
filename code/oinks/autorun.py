"""
autorun.py
this is a little helper script to run and then animate the results from a simulation
using oinks.OrbitalSystem()
"""

import numpy as np
import os
import nbody_kerr as oinks
import ThesisTools as farts

test_state = np.array([[[11,0],[0,0.17]],[[12,1],[0,0.30]],[[13,0],[0,0.30]]])
test_masses = np.array([5,10,15])

sim = oinks.OrbitalSystem(3,3,init_state = test_state, masses = test_masses,save_dir = "./output")
sim.TimeEvolve(500,comments="looking at if the mass array is dealt with properly.")
sim.movie()
moviepath = "./"+sim.dirname+"/movie/out.mp4"
os.system("open {}".format(moviepath))