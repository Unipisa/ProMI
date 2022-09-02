# Syntax to run genMutation.py
python3 genMutation.py -f [path of the fasta file of wt protein] -o [path to the folder contain the output] -m [path of the list of mutation file]

Example, open terminal in this directory run: <br>
python3 genMutation.py -o genMut/ -m list.txt -f testsingle.fasta <br>

## Syntax of the file contain the list of the mutation, each line: <br>
[new_aa] [start] [stop] [seq_id] <br>
<br>
For example: <br>
M 1 5 1 <br>
K 7 19 2 <br>
means: <br>
- from position 1 to position 5 in the first sequence (in case the fasta file contains more than 1 sequences) change the wt amino acid to M  <br>
- from position 7 to 19 in the second sequence change the wt amino acid to K <br>


<br> 

# Syntax to run predictBatch.py <br>
python3 predictBatch.py -f [path to folder contains fasta files] -o [path to output folder] -a [path to alphafold folder] -t [max template date - optional]
For example: <br>
python3 predictBatch.py -f /home/t.pham1/Sequences/genMut -o /home/t.pham1/Results/test -a /home/t.pham1/alphafold -t 2008-03-31 <br>


# Syntax to run alignStructure.py <br>
python3 alignStructure.py -f [Path of the folder contains subfolders, each subfolder is a prediction of alphafold for a mutation, the ranked_0.pdb file is inside the subfolders] -o [Output path] -c [Path to ChimeraX running file] -p [Path to the pdb model that used to align the input structures] -d [Optional, for printing log] <br>
For example: <br>
python3 alignStructure.py -f /Users/phamgiang/Documents/Study/Master/Thesis/Workspace/Validate/HIV1_PROTEASE_1/1GNO_MUTATIONS -o /Users/phamgiang/Documents/Study/Master/Thesis/Workspace/Validate/HIV1_PROTEASE_1/MD -c /Applications/ChimeraX-1.4-rc2022.05.29.app/Contents/MacOS/ChimeraX -p /Users/phamgiang/Documents/Study/Master/Thesis/Workspace/Validate/HIV1_PROTEASE_1/Raw/1gno.pdb -d  <br>
<br>

## Syntax for the directory structure <br>
The root directory of all the computation (which is pass as agrument -f) contains: <br>
- A subfolder named 'Share' which contains all the files needed for all the ligands and also the ligands.txt file and the config_pdb2gmx.txt file that has information about force field and water model. If there are more than one species of ligand, each species should be in a separated subfolder.  <br>
- A list of subfolders which will be the working directory for each MD simulation of each mutated protein. In those subfolder there is needed a .pdb of the pure protein. <br>
