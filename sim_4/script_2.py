no = 2

import hoomd
import hoomd.md

hoomd.context.initialize("");

# Build a lattice of molecules
uc = hoomd.lattice.unitcell(
  N=3,
  a1=[3, 0, 0],
  a2=[0, 3, 0],
  a3=[0, 0, 3],
  dimensions=3,
);

# Create snapshot
snapshot = uc.get_snapshot();

# Set particle positions
snapshot.particles.position[:] = [
  [0,0,0], [0.4, 0.4, 0], [0.8, 0.8, 0]
]

# Assign partcile types
snapshot.particles.types = ['A', 'B'];
snapshot.particles.typeid[0:2] = 0;
snapshot.particles.typeid[2:3] = 1;

# Set bonds between particles
snapshot.bonds.resize(2);
snapshot.bonds.types = ['polymer1', 'polymer2'];
snapshot.bonds.group[:] = [[0, 1], [1, 2]];
snapshot.bonds.typeid[:] = [0, 1]

# Replicate system to make into polymer system
# ERROR if you make this 5,5,5, then the nl can't be computed...
snapshot.replicate(10, 10, 10);
system = hoomd.init.read_snapshot(snapshot);

nl = hoomd.md.nlist.cell();
lj = hoomd.md.pair.lj(r_cut=1.5, nlist=nl);
lj.pair_coeff.set('A', 'A', epsilon=1.0, sigma=1.0);
lj.pair_coeff.set('B', 'B', epsilon=1.5, sigma=1.5);
lj.pair_coeff.set('A', 'B', epsilon=0.5, sigma=0.5);

# Define bonds
harmonic = hoomd.md.bond.harmonic();
harmonic.bond_coeff.set('polymer1', k=20.0, r0=0);
harmonic.bond_coeff.set('polymer2', k=5.0, r0=0);

hoomd.md.integrate.mode_standard(dt=0.01);
all = hoomd.group.all();
integrator = hoomd.md.integrate.nve(group=all);
integrator.randomize_velocities(kT=10.0, seed=42)

hoomd.analyze.log(
  filename=f"log-output_{no}.log",
  quantities=['potential_energy', 'temperature'],
  period=500,
  overwrite=True
);

hoomd.dump.gsd(f"trajectory_{no}.gsd", period=1, group=all, overwrite=True);
hoomd.run(100);
