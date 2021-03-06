# Example for: model.reorder_atoms()

# This will standardize the order of atoms in the model.

import sys

from modeller import *
from modeller.scripts import complete_pdb

env = environ()
env.io.atom_files_directory = ['.']
# Order the atoms according to a topology library:
env.libs.topology.read(file='/Library/modeller-9.9/modlib/top_heav.lib')
env.libs.parameters.read(file='/Library/modeller-9.9/modlib/par.lib')
mdl = model(env)
mdl.read(file=sys.argv[1])
mdl.reorder_atoms()
mdl.write(file=sys.argv[2])
