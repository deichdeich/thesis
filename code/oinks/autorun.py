import nbody_kerr as oinks
import ThesisTools as farts

#test_state = [[[7,0],[0,0.17]]]

sim = oinks.OrbitalSystem(300,1,interaction=None)
sim.TimeEvolve(3000,comments="no interaction")
#sim.trajplot()
sim.movie()