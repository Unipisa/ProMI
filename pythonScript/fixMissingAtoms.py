from modeller import *
from modeller.automodel import *
from modeller.scripts import complete_pdb
import os 

env = Environ()
env.io.atom_files_directory = ["../atom_files"]
log.verbose()
env.libs.topology.read(file="$(LIB)/top_heav.lib")
env.libs.parameters.read(file="$(LIB)/par.lib")

infolder = '/Users/phamgiang/Documents/Study/Master.nosync/Thesis/Workspace.nosync/scripts/exp_data/nofix2/3um8.pdb'
outfolder = '/Users/phamgiang/Documents/Study/Master.nosync/Thesis/Workspace.nosync/scripts/exp_data/fix/fix_3um8.pdb'

# files = os.scandir(infolder)
# for file in files:
#     if '.pdb' in file.name:
#         filename = os.path.join(infolder, file.name)
#         mdl = complete_pdb(env, filename)
#         temname = os.path.join(outfolder, 'tem_' + file.name)
#         outname = os.path.join()
#         mdl.write(temname)
#         infile = open(filename, 'r')
#         inlines = infile.readlines()

#         temfile = open(temname, 'r')
#         temlines = infile.readlines()
#         noend = temlines[0:len(temlines) - 1] # omit the last line "END"

#         for line in inlines:
#             if line[0:6] == 'HETATM':
#                 noend.append(line)
        
#         noend.append(temlines[-1])

#         outfile = open(os.path.join(outfolder, "fix_" + file.name), 'w')
#         outfile.write(''.join(noend))
#         outfile.close()



mdl = complete_pdb(env, infolder)
mdl.write(outfolder) 