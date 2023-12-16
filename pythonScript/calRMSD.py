from biopandas.pdb import PandasPdb
import os
import sys
import numpy as np
import pandas as pd 

# np.set_printoptions(threshold=sys.maxsize)

def smithwaterman (str1, str2, debug=False):
    nrow = len(str1)
    ncol = len(str2)

    mat = np.zeros([nrow+1, ncol+1], dtype=int) # len1 + 1 rows and len2 + 1 columns 
    for i in range(1,nrow+1):
        for j in range(1,ncol+1):
            if str1[i-1] == str2[j-1]:
                temlist = [max(0, mat[i-1,j] - 3), max(0, mat[i][j-1] - 3), max(0, mat[i-1][j-1] + 1)]
                mat[i][j] = max(temlist)
            else:
                temlist = [max(0, mat[i-1,j] - 3), max(0, mat[i][j-1] - 3), max(0, mat[i-1][j-1] - 1)]
                mat[i][j] = max(temlist)
    # print(mat)
    maxele = np.max(mat)
    maxidx = np.argmax(mat)
    rowmax = round(maxidx/(ncol+1)) - 1
    colmax = maxidx%(ncol+1) 

    
    i = rowmax
    j = colmax
    # pairs = [(i-1, j-1)]
    pairs1 = dict()
    pairs2 = dict() 
    
    while(i > 1 and j > 1):
        if mat[i-1][j-1] == max(mat[i-1][j-1], mat[i-1][j], mat[i][j-1]):
            i = i - 1
            j = j - 1
            pairs1[str(i)] = j
            pairs2[str(j)] = i
            # pairs.append((i-1, j-1))
        elif mat[i-1][j] == max(mat[i-1][j-1], mat[i-1][j], mat[i][j-1]):
            i = i - 1
            # pairs.append((i-1,j-1))
            pairs1[str(i)] = j
            pairs2[str(j)] = '-'
        else:
            j = j - 1
            # pairs.append((i-1, j-1))
            pairs1[str(i)] = '-'
            pairs2[str(j)] = i
        
        if mat[i][j] == 0:
            break
    return pairs1, pairs2


def getseq(atomdf, debug=False):
    codes = {"ALA": 'A', "ARG": "R", "ASN": 'N', "ASP": "D", 'ASX': 'B', 'CYS': 'C', 'GLU': 'E', \
     'GLN': 'Q', 'GLX': 'Z', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K', \
     'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', \
     'VAL': 'V', 'HIE': 'H', 'HID': 'H', 'HIP': 'H', 'ASH': 'D', 'GLH': 'E', 'CYM': 'C', \
     'CYX': 'C', 'LYN': 'K'}
    
    seen = set() 
    seq = ''

    maxlen = atomdf.tail(1)['residue_number']
    print(maxlen.tolist()[0])
    
    i = 0 # sometimes res_id doesn't start with 1, this variable is to pad '_' until the sequence really start
    for idx, atom in atomdf.iterrows():
        res_id = atom['residue_number']
        res_name = atom['residue_name']
        if res_id in seen:
            continue
        else: 
            seen.add (res_id)
            i += 1 
            while i < res_id: # add '_' to fill all the missing residues 
                seq += '_'
                i += 1
            seq += codes[res_name.upper()]
            
    # print(seen)
    return seq

            
def apply(path_exp, path_align, path_out, ligname, dist = 5, debug=True):
    if debug:
        print("Load file {}".format(path_align))
    pr_align = PandasPdb().read_pdb(path_align)

    # al_seq = ''.join(pr_align.amino3to1()['residue_name'].tolist())

    if debug:
        print("Load file {}".format(path_exp))
    pr_exp = PandasPdb().read_pdb(path_exp)


    # get all kind of ligands 
    ligands = pr_exp.df['HETATM']

    # get ligand of interest
    ligand = ligands[ligands['residue_name'] == ligname]
    # print(ligand)

    # list of residue id of alphafold structure near the ligand 
    al_resids = []

    # list of residue id of experimental structure near the ligand 
    ex_resids = []

    print(ligand.index)
    for lidx in ligand.index:
        # get atoms within radius = dist from ligand of experimental structure 
        l = ligands.iloc[lidx] 
        # print(l) 
        l_coor = (l.loc['x_coord'], l.loc['y_coord'], l.loc['z_coord'])

        al_dist = pr_align.distance(xyz = l_coor, records=('ATOM'))
        lig_atom_within = pr_align.df['ATOM'][al_dist < dist]
        al_resids += lig_atom_within['residue_number'].drop_duplicates().tolist()


        ex_dist = pr_exp.distance(xyz = l_coor, records=('ATOM'))
        ex_atom_within = pr_exp.df['ATOM'][ex_dist < dist]
        ex_resids += ex_atom_within['residue_number'].drop_duplicates().tolist()
    
    # list of residue id of atoms within the radius
    al_resids = list(set(al_resids))
    ex_resids = list(set(ex_resids))
    al_resids.sort()
    ex_resids.sort()
    print(al_resids)
    print(ex_resids)


    atom_pr_align = pr_align.df['ATOM']
    atom_pr_ex = pr_exp.df['ATOM']

    al_seq = getseq(atom_pr_align)
    ex_seq = getseq(atom_pr_ex)

    print(al_seq)
    print(ex_seq)

    pairs1, pairs2 = smithwaterman(al_seq, ex_seq)
    
    print(pairs1)
    print(pairs2)
    
    matchpairs = set()

    # take all the residues of al within radius and find the corresponding residue of ex
    for al_res in al_resids:
        try:
            cor_exres = pairs1[str(al_res)] # corresponding residue in experimental structure 
            matchpairs.add((al_res, cor_exres))

            # print(al_res, cor_exres)
            # print(al_seq[al_res], ex_seq[cor_exres])
            if debug:
                al_atoms1 = atom_pr_align[atom_pr_align['residue_number'] == al_res]
                print(al_res, al_atoms1['residue_name'].tolist())
                
                if cor_exres != '_':
                    ex_atoms1 = atom_pr_ex[atom_pr_ex['residue_number'] == cor_exres]
                    print(cor_exres, ex_atoms1['residue_name'].tolist())
                else:
                    print('_')
                print('\n')
        except:
            print("Not found {} in 1st aligned pairs".format(al_res))
    

    # take all the residues of ex within radius and find the corresponding residue of al 
    for ex_res in ex_resids:
        try:
            cor_alres = pairs2[str(ex_res)] # corresponding residue in alphafold structure 
            if (cor_alres, ex_res) not in matchpairs:
                matchpairs.add(cor_alres, ex_res)
            # print(ex_res, cor_alres)
            # print(al_seq[al_res], ex_seq[cor_exres])
            if debug:
                ex_atoms2 = atom_pr_ex[atom_pr_ex['residue_number'] == ex_res]
                print(ex_res, ex_atoms2['residue_name'].tolist())
                
                if cor_alres != '_':
                    al_atoms2 = atom_pr_align[atom_pr_align['residue_number'] == cor_alres]
                    print(cor_alres, al_atoms2['residue_name'].tolist())
                else:
                    print('_')
                print('\n')
        except:
            print("Not found {} in 2nd aligned pair".format(ex_res))
    
    if debug:
        print(matchpairs)

    # two lists to save al atoms and ex atoms to calculate RMSD 
    sc_al_atoms = []
    sc_ex_atoms = []

    bb_al_atoms = []
    bb_ex_atoms = []

    all_al_atoms = []
    all_ex_atoms = []

    backbone = ['CA','HA','N','C','O','HN','H']

    # for each matching pair of residue
    # take the atom and save to a new dataframe to be able to calculate RMSD 
    for pair in matchpairs:
        alresid, exresid = pair[0], pair[1] # get pair of matching residue ids

        al_atoms = atom_pr_align[atom_pr_align['residue_number'] == alresid]
        ex_atoms = atom_pr_ex[atom_pr_ex['residue_number'] == exresid]

        # print(al_atoms['residue_name'])
        # print(ex_atoms['residue_name'])

        # check correct match by checking name of residue, in case of mutation
        if al_atoms.shape[0] <= 0 or ex_atoms.shape[0] <= 0:
            continue
        else:
            if al_atoms['residue_name'].iloc[0] != ex_atoms['residue_name'].iloc[0]:
                print("Encounter {} and {}".format(al_atoms['residue_name'].iloc[0],\
                                                    ex_atoms['residue_name'].iloc[0]))

        for idx, ex_atom in ex_atoms.iterrows():
            if ex_atom['atom_name'] not in backbone:
                # get dataframe of 1 row of side chain
                corsc_al_atom = al_atoms[al_atoms['atom_name'] == ex_atom['atom_name']] 
                for idx2, al_atom1 in corsc_al_atom.iterrows():
                    sc_al_atoms.append(al_atom1)
                    sc_ex_atoms.append(ex_atom)
            else:
                # get dataframe of 1 row of backbone
                corbb_al_atom = al_atoms[al_atoms['atom_name'] == ex_atom['atom_name']]
                for idx3, al_atom2 in corbb_al_atom.iterrows():
                    bb_al_atoms.append(al_atom2)
                    bb_ex_atoms.append(ex_atom)
            
            # get dataframe of 1 row of all atoms excluding hydrogen 
            if len(ex_atom['atom_name']) > 0 and ex_atom['atom_name'][0] != 'H': # excluding hydrogen 
                cor_al_atom = al_atoms[al_atoms['atom_name'] == ex_atom['atom_name']]
                for idx4, al_atom3 in cor_al_atom.iterrows():
                    all_al_atoms.append(al_atom3)
                    all_ex_atoms.append(ex_atom)
            




    if debug:  
        print("Len of sidechain dataframe of al is {} and ex is {}".format(len(sc_al_atoms), \
                                                                           len(sc_ex_atoms)))
        print("Len of backbone dataframe of al is {} and ex is {}".format(len(bb_al_atoms), \
                                                                           len(bb_ex_atoms)))
        print("Len of all dataframe of al is {} and ex is {}".format(len(all_al_atoms), \
                                                                           len(all_ex_atoms)))
        print("\n")
        

    df_al_sc = pd.DataFrame(sc_al_atoms)
    df_ex_sc = pd.DataFrame(sc_ex_atoms)

    df_al_bb = pd.DataFrame(bb_al_atoms)
    df_ex_bb = pd.DataFrame(bb_ex_atoms)

    df_al_all = pd.DataFrame(all_al_atoms)
    df_ex_all = pd.DataFrame(all_ex_atoms) 

    rmsd_sc = PandasPdb.rmsd(df_al_sc, df_ex_sc, s=None)
    rmsd_bb = PandasPdb.rmsd(df_al_bb, df_ex_bb, s=None)
    rmsd_all = PandasPdb.rmsd(df_al_all, df_ex_all, s=None)

    print ("RMSD of sidechaine is {} over {} residues".format(rmsd_sc, len(matchpairs)))
    print ("RMSD of backbone is {} over {} residues".format(rmsd_bb, len(matchpairs)))
    print ("RMSD of all is {} over {} residues".format(rmsd_all, len(matchpairs)))

    # write output file 
    file_out = open(path_out, 'w+')
    text = 'AlphaFold2 file:   ' + path_align + '\n'
        
    text += 'Local RMSD:    ' + str(rmsd_all) + '\n'        
    text += 'Local RMSD sidechains: ' + str(rmsd_sc)+'\n'
    text += 'Local RMSD backbone: ' + str(rmsd_bb)+'\n\n'

    file_out.write(text)
    file_out.close()
    print("Save result file at {}".format(path_out))


    

    



    
    









