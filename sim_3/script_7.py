no = 7
# v5: Made all of them of type A
# v6
# tried making forces weaker so that they would be slower
# also removed bond from center so that we have double the number of polymers
# then made half of them A, the other half B
# v6
# Reduced forced even more; redcuced temperature by half; run for a longer time period

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

# Assign types of particles (0 is A, 1 is B)
snapshot.particles.typeid[0:5]=0;
snapshot.particles.typeid[5:10]=1;

# Create bonds between sets of paricles
snapshot.bonds.resize(8);
snapshot.bonds.group[:] = [
  [0, 1], [1, 2], [2, 3],
  [3, 4], [5, 6],
  [6, 7], [7, 8], [8, 9]
];

# Replicate the snapshot
snapshot.replicate(1, 10, 10);

# Initialize system
hoomd.init.read_snapshot(snapshot);

# Create neighbor list
nl = hoomd.md.nlist.cell();
dpd = hoomd.md.pair.dpd(r_cut=1.0, nlist=nl, kT=0.4, seed=1);
dpd.pair_coeff.set('A', 'A', A=1.0, gamma = 100.0);
dpd.pair_coeff.set('B', 'B', A=1.0, gamma = 100.0);
dpd.pair_coeff.set('A', 'B', A=1.0, gamma = 100.0);
nl.reset_exclusions(exclusions = []);

harmonic = hoomd.md.bond.harmonic();
harmonic.bond_coeff.set('polymer', k=5.0, r0=0);

hoomd.md.integrate.mode_standard(dt=0.01);
all = hoomd.group.all();
integrator = hoomd.md.integrate.nve(group=all);
integrator.randomize_velocities(kT=0.0, seed=42)

hoomd.analyze.log(
  filename=f"log-output_{no}.log",
  quantities=['potential_energy', 'temperature'],
  period=500,
  overwrite=True
);

hoomd.dump.gsd(f"trajectory_{no}.gsd", period=1, group=all, overwrite=True);
hoomd.run(1200);
