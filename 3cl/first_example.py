#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
# My first example with AutoDock Vina in python
#

from vina import Vina


v = Vina(sf_name='vina')

v.set_receptor('7rfs_receptor.pdbqt')

v.set_ligand_from_file('7rfs_ligand.pdbqt')
v.compute_vina_maps(center=[12.000, 0.000, 20.000], box_size=[60, 60, 60])
# print(1)
# Score the current pose
energy = v.score()
print('Score before minimization: %.3f (kcal/mol)' % energy[0])
print(1)

# Minimized locally the current pose
energy_minimized = v.optimize()
print('Score after minimization : %.3f (kcal/mol)' % energy_minimized[0])
v.write_pose('7rfs_ligand_minimized.pdbqt', overwrite=True)

# Dock the ligand
v.dock(exhaustiveness=32, n_poses=20)
v.write_poses('7rfs_ligand_vina_out.pdbqt', n_poses=5, overwrite=True)