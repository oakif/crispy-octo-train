# https://lost-contact.mit.edu/afs//umich.edu/user/j/o/joaander/Public/hoomd-web/doc/create_random_polymers-example.html

# ---- create_random_polymers.py ----
from hoomd import *
# from hoomd_script import *
import math
# parameters
phi_P = 0.25
n_poly = 600
T = 1.2
polymer1 = dict(bond_len=1.2, type=['A']*6 + ['B']*7 + ['A']*6,
                bond="linear", count=n_poly)
# perform some simple math to find the length of the box
N = len(polymer1['type']) * polymer1['count']
# generate the polymer system
init.create_random_polymers(box=data.boxdim(volume=math.pi * N / (6.0 * phi_P)), polymers=[polymer1],
                            separation=dict(A=0.35, B=0.35), seed=12)
# force field setup
harmonic = bond.harmonic()
harmonic.bond_coeff.set('polymer', k=330.0, r0=0.84)
lj = pair.lj(r_cut=3.0)
lj.pair_coeff.set('A', 'A', epsilon=1.0, sigma=1.0)
lj.pair_coeff.set('A', 'B', epsilon=1.0, sigma=1.0)
lj.pair_coeff.set('B', 'B', epsilon=1.0, sigma=1.0)
# dump every 100,000 steps
dump.xml(filename="create_random_polymers.xml", vis=True)
dump.dcd(filename="create_random_polymers.dcd", period=100000)
# integrate NVT for a bunch of time steps
all = group.all()
integrate.mode_standard(dt=0.005)
integrate.nvt(group=all, T=1.2, tau=0.5)
run(2000)
# uncomment the next run() command if you have a few hours to spare
# running this on a GPU the resulting dump files should show the
# polymers self-assembling into the hex phase
# run(10e6)
