import nbody_kerr as oinks
import ThesisTools as farts
import multiprocessing as mp

nparticles = 15500

init_vecs = farts.make_initial_conditions(nparticles)

def a1(t):
    if t < 15000:
        return 0
    else:
        return .99

a_list = [a1]

def a_chooser(inter):
    sim = oinks.OrbitalSystem(nparticles,3,init_state = init_vecs,a=inter,save_dir = "/Volumes/DRIVEYWIVEY/thesis/output")
    sim.TimeEvolve(40000, comments = "")
    sim.movie()
    

def mp_write():
    num_proc = 2
    pool = mp.Pool(num_proc)
    pool.map(a_chooser,a_list)
    