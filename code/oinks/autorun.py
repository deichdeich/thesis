
"""
autorun.py
this is a little helper script to run and then animate the results from a simulation
using oinks.OrbitalSystem()
"""


import nbody_kerr as oinks
import ThesisTools as farts

#test_state = [[[7,0],[0,0.17]]]

sim = oinks.OrbitalSystem(2,1)
sim.TimeEvolve(5000,comments="")
sim.movie()