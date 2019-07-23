no = 1
# v1 updated with lj interactions between molecules

import hoomd
import hoomd.md

hoomd.context.initialize("");

# Create box containing polymer with N particles
snapshot = hoomd.data.make_snapshot(
  N=10,
  box=hoomd.data.boxdim(Lx=10, Ly=0.5, Lz=0.5),
  particle_types=['A', 'B'],
  bond_types=['polymer']
);

# Assign positions of polymer
snapshot.particles.position[:] = [
  [-4.5, 0, 0], [-3.5, 0, 0],
  [-2.5, 0, 0], [-1.5, 0, 0],
  [-0.5, 0, 0], [0.5, 0, 0],
  [1.5, 0, 0], [2.5, 0, 0],
  [3.5, 0, 0], [4.5, 0, 0]
];

# Assign types to particles (0 is A, 1 is B)
snapshot.particles.typeid[0:5]=0;
snapshot.particles.typeid[5:10]=1;

# Create bonds between sets of paricles
snapshot.bonds.resize(8);
snapshot.bonds.group[:] = [
  [0, 1], [1, 2], [2, 3], [3, 4],
  [5, 6], [6, 7], [7, 8], [8, 9]
];

# Replicate the snapshot
snapshot.replicate(1, 20, 20);

# Initialize system
hoomd.init.read_snapshot(snapshot);

# Specify particle interactions
nl = hoomd.md.nlist.cell();
lj = hoomd.md.pair.lj(r_cut=1.0, nlist=nl);
lj.pair_coeff.set('A', 'A', epsilon=1.0, sigma=1.0);
lj.pair_coeff.set('B', 'B', epsilon=1.5, sigma=1.5);
lj.pair_coeff.set('A', 'B', epsilon=0.5, sigma=0.5);

# Specify bonded interactions
harmonic = hoomd.md.bond.harmonic();
harmonic.bond_coeff.set('polymer', k=5.0, r0=0);

# Integrator
hoomd.md.integrate.mode_standard(dt=0.01);
all = hoomd.group.all();
integrator = hoomd.md.integrate.nve(group=all);
integrator.randomize_velocities(kT=0.0, seed=42)

# Output

hoomd.analyze.log(
  filename=f"log-output_{no}.log",
  quantities=['potential_energy', 'temperature'],
  period=500,
  overwrite=True
);

hoomd.dump.gsd(f"trajectory_{no}.gsd", period=100, group=all, overwrite=True);
hoomd.run(10e5);