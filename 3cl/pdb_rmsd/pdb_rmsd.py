# import itertools
# import time
import pymol
pdb_list = ['6BHT']
for  i in pdb_list:
    pymol.cmd.fetch(i)
a = pymol.cmd.get_object_list()
# combinations = list(itertools.combinations(pdb_list, 2))
# print(combinations)

