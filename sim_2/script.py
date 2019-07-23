import hoomd
import hoomd.md

# Define polymer
polymer1 = dict(
  bond_len = 1.0,
  type = ['A'] * 32, # number of beads
  bond = 'linear',
  count = 1, # number of polymers
  seed = 5
)

# Find length of bond required; equal to the number of total monomers
N = len(polymer1['type']) * polymer1['count']

# Generate a system consisting of a single polymer in a random polymer configuration
system = hoomd.deprecated.init.create_random_polymers(
  box = hoomd.BoxDim(100, 100, 100), # large box
  polymers = [polymer1],
  separation = dict(A = 0.25, B = 0.25)
)

# Set position of first bead (particle) to the origin of box
system.particles[0].position = (0, 0, 0)

# Set up force field. It is essentially a harmonic spring hwich holds beads together
harmonic = hoomd.bond.harmonic()
harmonic.bond_coeff.set('polymer', k = 400.0, r0 = 1.0)

# Define pair force between monomers which attracts them to each other, but also ensures beads don't overlap
lj = hoomd.pair.lj(r_cut = 2**(1./6.))
lj.pair_coeff.set('A', 'A', epsilon = 1.0, sigma = 1.0)

# NVT integration
all = hoomd.group.all()
tail = hoomd.group.tags(tag_min = 1, tag_max = len(all) - tag_min) # keep beads undr tag_min fixed; integrate all beads from tag_min to end
hoomd.integrate.mode_standard(dt = 0.005)
hoomd.integrate.nvt(group = tail, T = 1.0, tau = 1.0) # tau is relaxation time?

hoomd.analyze.log(filename="log-output.log", quantities=['potential_energy'], period=100, overwrite=True) # log the stuff you want
hoomd.dump.gsd('trajectory.gsd', period=5, group=all, overwrite=True) # output trajectory

hoomd.run(1e4)