import nbody_kerr as oinks
import ThesisTools as farts

#test_state = [[[7,0],[0,0.17]]]

sim = oinks.OrbitalSystem(2,1)
sim.TimeEvolve(5000,comments="")
sim.movie()