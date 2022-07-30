# Syntax to run genMutation.py
python3 genMutation.py -f [path of the fasta file of wt protein] -o [path to the folder contain the output] -m [path of the list of mutation file]

Example, open terminal in this directory run:
python3 genMutation.py -o genMut/ -m list.txt -f testsingle.fasta 

## Syntax of the file contain the list of the mutation, each line:
[new_aa] [start] [stop] [seq_id]

For example: 
M 1 5 1
K 7 19 2
means:
- from position 1 to position 5 in the first sequence (in case the fasta file contains more than 1 sequences) change the wt amino acid to M 
- from position 7 to 19 in the second sequence change the wt amino acid to K




# Syntax to run predictBatch.py 
python3 predictBatch.py -f [path to folder contains fasta files] -o [path to output folder] -a [path to alphafold folder] -t [max template date - optional]
For example: 
python3 predictBatch.py -f /home/t.pham1/Sequences/genMut -o /home/t.pham1/Results/test -a /home/t.pham1/alphafold -t 2008-03-31


# Syntax to run alignStructure.py
python3 alignStructure.py -f [Path of the folder contains subfolders, each subfolder is a prediction of alphafold for a mutation, the ranked_0.pdb file is inside the subfolders] -o [Output path] -c [Path to ChimeraX running file] -p [Path to the pdb model that used to align the input structures] -d [Optional, for printing log]
For example:
python3 alignStructure.py -f /Users/phamgiang/Documents/Study/Master/Thesis/Workspace/Validate/HIV1_PROTEASE_1/1GNO_MUTATIONS -o /Users/phamgiang/Documents/Study/Master/Thesis/Workspace/Validate/HIV1_PROTEASE_1/MD -c /Applications/ChimeraX-1.4-rc2022.05.29.app/Contents/MacOS/ChimeraX -p /Users/phamgiang/Documents/Study/Master/Thesis/Workspace/Validate/HIV1_PROTEASE_1/Raw/1gno.pdb -d 

# Syntax to run formComplexBatch.py 
python3 formComplexBatch.py -f [Path the the folder contains subfolders, each subfolder is for one protein] -c [Path to the config file which has information of ligands]

## Syntax for the directory structure
The root directory of all the computation (which is pass as agrument -f) contains:
- A subfolder named 'Share' which contains all the files needed for all the ligands and also the ligands.txt file and the config_pdb2gmx.txt file that has information about force field and water model. If there are more than one species of ligand, each species should be in a separated subfolder. 
- A list of subfolders which will be the working directory for each MD simulation of each mutated protein. In those subfolder there is needed a .pdb of the pure protein.


## Syntax of the ligands.txt file 
[number_types_ligand]
[======] # separated line
[name_of_the_ligand]
[number_ligands_of_type_1] 
[Path_to_list_of_.gro_files_of_ligand_type_1]
[Path_to_.itp_file_of_ligand_type_1]
[Path_to_.prm_file_of_ligand_type_1]
[Path_to_posre_file_of_ligand_type_1]
[======]
[name_of_the_ligand]
[number_ligands_of_type_2] 
[Path_to_list_of_.gro_files_of_ligand_type_2]
[Path_to_.itp_file_of_ligand_type_2]
[Path_to_.prm_file_of_ligand_type_2]
[Path_to_posre_file_of_ligand_type_2]
....
[======]
[name_of_the_ligand]
[number_ligands_of_type_n] 
[Path_to_list_of_.gro_files_of_ligand_type_n]
[Path_to_.itp_file_of_ligand_type_n]
[Path_to_.prm_file_of_ligand_type_n]
[Path_to_posre_file_of_ligand_type_n]

All the paths need to be absolute 
The name of the ligand needs to be the same with the name declared in the section [ moleculetype ] of the file .itp 
The separate lines are mandatory

## Syntax of the config.txt file 
pdb2gmx: [option_to_run_pdb2gmx]
editconf: [option_to_run_editconf]
genion: [option_to_run_genion]