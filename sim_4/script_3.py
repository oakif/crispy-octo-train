# File stuff
trajectory_number = 3

# Bond stuff
polymer_A_beads = 8
polymer_B_beads = 8
bond_length = 1.25
polymer_separation = 2
polymer_length_dimension = 5
polymer_side_dimension = 10

# Calculated variables
total_beads = polymer_A_beads + polymer_B_beads
polymer_A_length = polymer_A_beads * bond_length
polymer_B_length = polymer_B_beads * bond_length
max_polymer_length = max(polymer_A_length, polymer_B_length)
total_bonds = polymer_A_beads + polymer_B_beads - 2

#######

import hoomd
import hoomd.md

hoomd.context.initialize("");

# Build a lattice of molecules
uc = hoomd.lattice.unitcell(
  N = total_beads,
  a1 = [max_polymer_length + polymer_separation, 0, 0],
  a2 = [0, 2 * polymer_separation, 0],
  a3 = [0, 0, 2 * polymer_separation],
  dimensions = 3,
);

# Create snapshot
snapshot = uc.get_snapshot();

# Set particle positions

polymer_A_positions = [[i * bond_length, - polymer_separation / 2, 0] for i in range(polymer_A_beads)]
polymer_B_positions = [[i * bond_length, polymer_separation / 2, 0] for i in range(polymer_B_beads)]

position_list = polymer_A_positions + polymer_B_positions
for position in position_list: position[0] -= max_polymer_length / 2

snapshot.particles.position[:] = position_list

# Assign partcile types

snapshot.particles.types = ['A', 'B'];

snapshot.particles.typeid[0:polymer_A_beads] = 0;
snapshot.particles.typeid[polymer_A_beads:total_beads] = 1;

# Set bonds between particles

snapshot.bonds.resize(total_bonds);
snapshot.bonds.types = ['polymer1', 'polymer2'];

polymer_A_bonds = [[i, i + 1] for i in range(polymer_A_beads - 1)]
polymer_A_bond_types = [0 for i in range(polymer_A_beads - 1)]
polymer_B_bonds = [[polymer_A_beads + i, polymer_A_beads + i + 1] for i in range(polymer_B_beads - 1)]
polymer_B_bond_types = [1 for i in range(polymer_B_beads - 1)]

polymer_bonds = polymer_A_bonds + polymer_B_bonds
polymer_bond_types = polymer_A_bond_types + polymer_B_bond_types

snapshot.bonds.group[:] = polymer_bonds
snapshot.bonds.typeid[:] = polymer_bond_types

# Replicate system to make into polymer system

snapshot.replicate(polymer_length_dimension, polymer_side_dimension, polymer_side_dimension);
system = hoomd.init.read_snapshot(snapshot);

# Assign LJ potentials

nl = hoomd.md.nlist.cell();
lj = hoomd.md.pair.lj(r_cut=1.5, nlist=nl);
lj.pair_coeff.set('A', 'A', epsilon=1.0, sigma=1.0);
lj.pair_coeff.set('B', 'B', epsilon=1.0, sigma=1.0);
lj.pair_coeff.set('A', 'B', epsilon=1.0, sigma=1.0);

# Define bonds
harmonic = hoomd.md.bond.harmonic();
harmonic.bond_coeff.set('polymer1', k=250.0, r0=0);
harmonic.bond_coeff.set('polymer2', k=400.0, r0=0);

hoomd.md.integrate.mode_standard(dt=0.001);
all = hoomd.group.all();
integrator = hoomd.md.integrate.nve(group=all);
integrator.randomize_velocities(kT=250.0, seed=42)

hoomd.analyze.log(
  filename=f"log-output_{trajectory_number}.log",
  quantities=['potential_energy', 'temperature'],
  period=1,
  overwrite=True
);

hoomd.dump.gsd(f"trajectory_{trajectory_number}.gsd", period=1, group=all, overwrite=True);
hoomd.run(1e3);