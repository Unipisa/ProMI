import argparse as ap
from curses import meta
from genericpath import isfile
import json
import os
import shutil


def parseAgruments():
    parser = ap.ArgumentParser(description='Generate mutated sequences') 
    parser.add_argument('-f', '--meta', type=str, nargs='?', \
        help='Path to meta file', required=True)
    parser.add_argument('-p', '--pdb', type=str, nargs='?', \
        help='Path to pdb file of wild type protein with ligand(s)', required=True)
    parser.add_argument('-o', '--outpath', type=str, nargs='?', \
        help='Path of output directory which will contain the "Share" folder', required=True)
    parser.add_argument('-i', '--data', type=str, nargs='?', \
        help='Path to data folder', required=True)
    parser.add_argument('-c', '--chimera', type=str, nargs='?', \
        help='Path to chimeraX running file', default='/user/bin/chimerax')
    parser.add_argument('-d', '--debug', \
        help="Print log to debug or not", action='store_true')
    return parser

# read meta data and save config file from the provided information 
def parseMetaFile(metapath, outpath, debug=False):
    meta = {'type': 'homo', 'ff': 'amber', 'num_chain': 1, 'lig': '', \
        'ff_name': "amber99sb", 'water': 'tip3p', 'bt': "dodecahedron", \
        'p_name': "NA", 'n_name': "CL"} # default meta data 
    if not os.path.isfile(metapath):
        print("The path of the pdb file is not correct")
        exit()
    file = open(metapath, 'r')
    lines = file.readlines() 
    for line in lines:
        line = line.strip() 
        if len(line) >= 1 and line[0] != '#': # omit command line marked by '#'
            words = line.split('=')
            meta[words[0].strip()] = words[1].strip()
    try:
        meta['num_chain'] = int(meta['num_chain'])
    except:
        print("The number of chain in meta file is not a number")
        exit()

    if meta['type'] == 'mono':
        meta['num_chain'] = 1
    if debug: 
        print(meta)
    # generate config file from meta file 
    config = 'pdb2gmx: -ff ' + meta['ff_name'] + ' -water ' + meta['water'] + ' -ignh\n'
    config += 'edifconf: -bt ' + meta['bt']+ ' -d 1.0\n'
    config += 'genion: -pname ' + meta['p_name'] + ' -nname ' + meta['n_name'] + '\n'
    config += 'grompp: -maxwarn 50\n'

    if not os.path.isdir(os.path.join(outpath, 'Share')):
        if debug:
            print("Creating 'Share' directory")
        os.mkdir(os.path.join(outpath, 'Share'))
    fname = os.path.join(outpath, 'Share','config.txt')
    fconfig = open(fname, 'w')
    fconfig.write(config)
    fconfig.close()

    return meta

def parseLine(line, debug=False):
    key = str(line[0:6]).strip()
    id = int(''.join(line[6:11]))
    atom = str(line[11:17]).strip()
    name = str(line[17:21]).strip()
    chain = str(line[21:22]).strip()
    resid = int(line[22:26])
    if key != 'TER':
        x = float(line[26:38])
        y = float(line[38:46])
        z = float(line[46:54])
        a = float(line[54:60])
        b = float(line[60:66])
        kind = str(line[66:len(line)])
    else:
        x, y, z, a, b = 0, 0, 0, 0, 0
        kind = ''

    if debug:
        print(key, id, atom, name, chain, resid, x, y, z, a, b, kind)
    linedict = {'key': key, \
        'aid': id, 'atom': atom ,'name': name, \
        'chain': chain, 'rid': resid, \
        'line': line}
    return linedict


def getLigands(pdbpath, debug=False):
    if not os.path.isfile(pdbpath):
        print("The path of the pdb file is not correct")
        exit()
    file = open(pdbpath, 'r')
    lines = file.readlines() 

    # dictionary of list: {'chain': [list_of_lines]}
    prochains = {} 
    # dictionary (key is chain) of dictionary (key is resid) of dictionary (key is name of ligand) of list:
    # {'chainA': {'resid1': {'name': name, 'lines': [lines]}, 'resid2': {'name': name, 'lines': [lines]}}, \
    #  'chainB': {'resid1': {'name': name, 'lines': [lines]}, 'resid2': {'name': name, 'lines': [lines]}}
    ligs = {} 

    for line in lines:
        if line[0:4] == 'ATOM' or line[0:3] == 'TER':
            # parse ATOM line
            linedict = parseLine(line)
            chain = linedict['chain']
            if chain not in prochains:
                if debug:
                    print("New protein chain: {}".format(chain))
                prochains[chain] = [linedict['line']]
            else:
                # if debug:
                #     print("Append line to old chain {}".format(chain))
                prochains[chain].append(linedict['line'])
        if line[0:6] == 'HETATM' and 'HOH' not in line:
            # parse HETATM line
            linedict = parseLine(line)
            chain = linedict['chain']
            rid = linedict['rid'] # get residue if of ligand
            name = linedict['name'] # get the name of the ligand
            if chain not in ligs:
                if debug:
                    print("New ligand chain {}".format(chain))
                ligs[chain] = {}
                if debug: 
                    print("New residue with id {} for ligand of chain {}".format(rid, chain))
                ligs[chain][rid] = {'name': name, 'lines': [linedict['line']]}
            else: # in case the same old chain
                # if debug:
                    # print("Old ligand chain {}".format(chain))
                if rid not in ligs[chain]:
                    if debug: 
                        print("New ligand {} of chain {}".format(rid, chain))
                    ligs[chain][rid] = {'name': name, 'lines': [linedict['line']]}
                else:
                    # if debug:
                        # print("Add new line to old ligand {} of old chain {}".format(rid, chain))
                    ligs[chain][rid]['lines'].append(linedict['line'])
    return prochains, ligs

# Function to read the metadata to see how many ligands and ions to consider
# Then create the ligands.txt file
# Create one file for each ligand and then add hydrogen by calling chimera     
# Also prepare coupling group for .mdp files 
def separateLigands(prodict, ligdict, metadict, ionsdict, outpath, datapath, debug=False):
    # get the infor from the metadata dictionary 
    num_chain = metadict['num_chain']
    consider = metadict['lig']
    protype = metadict['type']

    if num_chain > len(prodict.keys()):
        print("The number of chain in the meta file excesses the number of chain in the pdb file")
        print("Check again")
        exit()
    if num_chain < len(prodict.keys()):
        if protype == 'homo':
            print("Will consider only {} chain(s) in total of {} chains of a homogenous protein".format(num_chain, len(prodict.keys())))
        if protype == 'complex':
            print("The protein is not homogenous \
                but the number of chain to consider is less than the number of chain from the pdb file \
                please check again")
            exit()

    if not os.path.isdir(outpath):
        print("Outpath {} is not a directory".format(outpath))
        exit()

    if not os.path.isdir(os.path.join(outpath, 'Share')):
        if debug:
            print("Creating 'Share' directory")
        os.mkdir(os.path.join(outpath, 'Share'))

    # general information (name, number of ligands) about all the ligands
    liginfor = {'ligands':{}, 'ions': {}} # {'ligands': {name: [number of that type]}, 'ions':{name : [number of that type]}}

    ligscons = {} # list of ligands of consider {'nameid': [line]}
    ionscons = {} # list of ions of consider {'nameid': [line]}
    savedlig = [] # list of filename of saved ligands to return 
    for count, chain in enumerate(ligdict.keys()): 
        count += 1
        if count > num_chain:
            break # only consider number of chain enter from the meta file for homogenous protein
        for id in ligdict[chain].keys():
            # each chain and id give a unique ligand with name and lines
            name = ligdict[chain][id]['name'] 
            num_line = len(ligdict[chain][id]['lines'])

            if num_line > 1: # not an ion
                if name not in liginfor['ligands']: 
                    liginfor['ligands'][name] = 1
                    nameid = name
                    ligscons[nameid] = ligdict[chain][id]['lines'] 
                    fname = os.path.join(outpath, 'Share', nameid + '.pdb')
                else: 
                    nameid = name + '_' + str(liginfor['ligands'][name])
                    ligscons[nameid] = ligdict[chain][id]['lines']
                    fname = os.path.join(outpath, 'Share', nameid + '.pdb')
                    liginfor['ligands'][name] += 1 # this must be after nameid 
                if debug:
                    print("Save ligand to file {}".format(fname))
                f = open(fname, 'w')
                f.write(''.join(ligdict[chain][id]['lines']))
                f.close()
                savedlig.append(nameid + '.pdb')
            else: # an ion with only one atom 
                if name in ionsdict: # only consider ion supported by the force field 
                    newlines = ligdict[chain][id]['lines'] # this is a list of string

                    # change the name of the ion if the name is not like the name declared in the force field
                    if name != ionsdict[name]:
                        if debug:
                            print("Change the name of ion {} to {}".format(name, ionsdict[name]))
                        newname = ''
                        for i in range (4 - len(ionsdict[name])):
                            newname += ' '
                        newname += ionsdict[name]
                        # find old name to replace
                        oldname = ''
                        for i in range(4 - len(name)):
                            oldname += ' '
                        oldname += name

                        for i, subline in enumerate(newlines):
                            subline = subline.replace(oldname, newname, 1)
                            newlines[i] = subline

                    if name not in liginfor['ions']:
                        liginfor['ions'][name] = 1
                        nameid = name
                        ionscons[nameid] = newlines
                        # fname = os.path.join(outpath, 'Share', nameid + '.pdb')
                    else:
                        nameid = name + '_' + str(liginfor['ions'][name]) 
                        ionscons[nameid] = newlines
                        # fname = os.path.join(outpath, 'Share', nameid + '.pdb')
                        liginfor['ions'][name] += 1 # this must be after nameid
                # there is not else, omit all kind of ions that is not supported by the forcefield 
    if debug:
        print(liginfor)
    
    if consider not in liginfor['ligands']:
        print("The ligand of interest {} is not in the ligand list".format(consider))
        exit()
    
    # creating ligands.txt file 
    if debug:
        print("Create ligands.txt file")
        
    ligfilelines = [str(len(liginfor['ligands'].keys())) + '\n'] # first line is the number of ligands
    ligfilelines.append('=============\n') # separate line
    ligfilelines.append(consider + '\n') # name of ligand of interest 
    # each ligand of that type add one .gro file
    ligfilelines.append(str(liginfor['ligands'][consider])+ '\n') # number of ligand belongs to species declared above 
    for i in range (liginfor['ligands'][consider]):
        if i == 0:
            line = '../Share/' + consider + '.acpype/' + consider + '_fix.gro\n'
        else:
            line = '../Share/' + consider + '_' + str(i) + '.acpype/' + consider + '_' + str(i) + '_fix.gro\n'
        ligfilelines.append(line)
    ligfilelines.append('../Share/' + consider + '.acpype/' + consider + '_fix.itp\n')
    ligfilelines.append('../Share/' + consider + '.acpype/' + consider + '_fix.prm\n')
    ligfilelines.append('../Share/' + consider + '.acpype/' + 'posre_' + consider + '.itp\n')
    # add other ligand 
    for lig in liginfor['ligands'].keys():
        if lig != consider: # already add considered ligand above
            ligfilelines.append('=============\n') # separate line
            ligfilelines.append(lig + '\n') # name of ligand 
            ligfilelines.append(str(liginfor['ligands'][lig])+ '\n')
            for i in range (liginfor['ligands'][lig]):
                if i == 0:
                    line = '../Share/' + lig + '.acpype/' + lig + '_fix.gro\n'
                else:
                    line = '../Share/' + lig + '_' + str(i) + '.acpype/' + lig + '_' + str(i) + '_fix.gro\n'
                ligfilelines.append(line)
            ligfilelines.append('../Share/' + lig + '.acpype/' + lig + '_fix.itp\n')
            ligfilelines.append('../Share/' + lig + '.acpype/' + lig + '_fix.prm\n')
            ligfilelines.append('../Share/' + lig + '.acpype/' + 'posre_' + lig + '.itp\n')
    ligfile = open(os.path.join(outpath, 'Share', 'ligands.txt'), 'w')
    ligfile.write(''.join(ligfilelines))
    ligfile.close()

    
    # write ions.txt file if there is any ion(s)
    # the ions.txt file contain:
    # first line is the number of species of ions will be add 
    # second line and so on will be the line of each ion
    if liginfor['ions']:
        if debug:
            print("Found some ions, creating ions.txt file")
        ionsfileline = [str(len(liginfor['ions'].keys())) + '\n']
        for ionid in ionscons:
            ionsfileline.extend(ionscons[ionid])
            
        ionfile = open(os.path.join(outpath, 'Share', 'ions.txt'), 'w')
        ionfile.write(''.join(ionsfileline))
        ionfile.close()
    if debug:
        print(savedlig)


    # from here prepare coupling group of .mdp files 
    # get information about ions and ligands from liginfor 
    # only couple protein with ligand, not ion 
    group1 = 'Protein'
    if liginfor['ions']: # if there are some ions to be considered 
        for ionname in liginfor['ions'].keys():
            group1 += '_' + ionname
    if liginfor['ligands']: # if there are some ligands to be considered 
        # consider the ligand of interest first 
        group1 += '_' + consider
        for ligname in liginfor['ligands'].keys():
            if ligname != consider:
                group1 += '_' + ligname
    if debug:
        print(group1)

    group2 = 'non-Protein'
    if liginfor['ions']: # if there are some ions to be considered 
        for ionname in liginfor['ions'].keys():
            group2 += '_&_!' + ionname
    if liginfor['ligands']: # if there are some ligands to be considered 
        # consider the ligand of interest first 
        group2 += '_&_!' + consider
        for ligname in liginfor['ligands'].keys():
            if ligname != consider:
                group2 += '_&_!' + ligname
    
    if debug:
        print(group2)

    

    # now add group to template .mdp file from folder data 
    # file nvt.mdp, copy file from data folder then put to Share folder and modify coupling group
    if os.path.isfile(os.path.join(datapath, 'nvt.mdp')):
        ftnvt = open(os.path.join(datapath, 'nvt.mdp'))
        nvtlines = ftnvt.readlines()
        ftnvt.close()
        for nvtid, nvtline in enumerate(nvtlines):
            if 'tc-grps' in nvtline:
                nvtline = nvtline.replace('Protein_Ligands', group1)
                nvtlines[nvtid] = nvtline.replace('Water_and_ions', group2)
        fnvt = open(os.path.join(outpath, 'Share', 'nvt.mdp'), 'w')
        fnvt.write(''.join(nvtlines))
        fnvt.close()
    else:
        print("Cannot file nvt.mdp file in the data folder")

    # file npt.mdp 
    if os.path.isfile(os.path.join(datapath, 'npt.mdp')):
        ftnpt = open(os.path.join(datapath, 'npt.mdp'))
        nptlines = ftnpt.readlines()
        ftnpt.close()
        for nptid, nptline in enumerate(nptlines):
            if 'tc-grps' in nptline:
                nptline = nptline.replace('Protein_Ligands', group1)
                nptlines[nptid] = nptline.replace('Water_and_ions', group2)
        fnpt = open(os.path.join(outpath, 'Share', 'npt.mdp'), 'w')
        fnpt.write(''.join(nptlines))
        fnpt.close()
    else:
        print("Cannot file nvt.mdp file in the data folder")

    # file md.mdp 
    if os.path.isfile(os.path.join(datapath, 'md.mdp')):
        ftmd = open(os.path.join(datapath, 'md.mdp'))
        mdlines = ftmd.readlines()
        ftmd.close()
        for mdid, mdline in enumerate(mdlines):
            if 'tc-grps' in mdline:
                mdline = mdline.replace('Protein_Ligands', group1)
                mdlines[mdid] = mdline.replace('Water_and_ions', group2)
        fmd = open(os.path.join(outpath, 'Share', 'md.mdp'), 'w')
        fmd.write(''.join(mdlines))
        fmd.close()
    else:
        print("Cannot file npt.mdp file in the data folder")

    # copy em.mdp and ions.mdp from data to Share directory
    emtar = os.path.join(datapath, 'em.mdp')
    emdes = os.path.join(outpath, 'Share', 'em.mdp')
    if os.path.isfile(emtar):
        shutil.copy(emtar, emdes)
    else:
        print("Cannot find em.mdp in data folder")
    
    ionstar = os.path.join(datapath, 'ions.mdp')
    ionsdes = os.path.join(outpath, 'Share', 'ions.mdp')
    if os.path.isfile(ionstar):
        shutil.copy(ionstar, ionsdes)
    else:
        print("Cannot find ions.mdp in data folder")

    # copy mmpbsa.in from data to Share directory
    mmtar = os.path.join(datapath, 'mmpbsa.in')
    mmdes = os.path.join(outpath, 'Share', 'mmpbsa.in')
    if os.path.isfile(mmtar):
        shutil.copy(mmtar, mmdes)
    else:
        print("Cannot find mmpbsa.in in data folder")

    return savedlig

# Run ChimeraX container to add Hydrogen to ligands 
def processLigand(savedlig, outpath, chipath, debug=False):
    # All the file need to be process by chimera and acpype will be in "Share" folder
    # get current directory to be able to change back 
    curdir = os.getcwd()
    # Get absolute path of "Share" folder 
    shareabs = os.path.abspath (os.path.join(outpath, "Share"))

    acpypefile = '#!/bin/sh\n'
    acpypefile += 'WS=/tmp/workspace/\n'
    acpypefile += 'cd $WS\n'

    os.chdir("/") # cd to the root folder in case of macos
    for lig in savedlig:
        if debug:
            print("Adding hydrogen to ligand {}".format(lig))
        ligpath = shareabs + '/' + lig
        cxcfile = 'addh\n'
        cxcfile += 'save ' + ligpath + " models #1\n"
        cxcfile += 'exit\n'
        # if debug:
        #     print(cxcfile)
        cxcpath = os.path.join(shareabs, 'add_hydrogen.cxc')
        f = open(cxcpath, 'w')
        f.write(cxcfile)
        f.close()
        # run ChimeraX 
        command = ['.' + chipath, '--nogui', '--silent' ,ligpath, cxcpath]
        # if debug:
        #     print(' '.join(command))
        os.system(' '.join(command))
        acpypefile += 'acpype -f -i ' + lig + ' -o gmx -d\n'

        # delete .cxc file 
        # if debug:
        #     print("Remove {}".format(cxcpath))
        os.remove(cxcpath)

    facname = os.path.join(shareabs, "run_acpype.sh")
    fac = open(facname, 'w')
    fac.write(acpypefile)
    fac.close()
    
    os.chdir(curdir)
    if debug:
        print("The current directory is {}".format(curdir))

def genInputGmx(outpath, debug=False):
    ligfname = os.path.join(outpath, 'Share', 'ligands.txt')
    ionname = os.path.join(outpath, "Share", 'ions.txt')
    
    nlig = 0
    if os.path.isfile(ligfname):
        flig = open(ligfname, 'r')
        try:
            nlig += int(flig.readline().strip())
        except:
            print("Cannot convert {} to number of species of ligands".format(flig.readline().strip()))
    
    nion = 0
    if os.path.isfile(ionname):
        fion = open(ionname, 'r')
        try:
            nion += int(fion.readline().strip())
        except:
            print("Cannot convert {} to number of species of ions".format(fion.readline().strip()))

    if debug:
        print("Number species of ions is {}".format(nion))
        print("Number species of ligands is {}".format(nlig))

    # create "Share" directory if not exist
    if not os.path.isdir(os.path.join(outpath, 'Share')):
        os.mkdir(os.path.join(outpath, 'Share'))
    
    # generate input_genion.txt 
    fgenionname = os.path.join(outpath, 'Share', 'input_genion.txt')
    fgenion = open(fgenionname, 'w')
    fgenion.write('SOL\n')
    fgenion.close()

    # generate input_coupling.txt
    fcouplename = os.path.join(outpath, 'Share', 'input_coupling.txt')
    fcouple = open(fcouplename, 'w')
    # the first group
    fcouplecont = '1'
    for i in range(nlig + nion):
        fcouplecont += ' | ' + str(13 + i)
    fcouplecont += '\n'
    # the second group
    fcouplecont += '11'
    for i in range(nlig + nion):
        fcouplecont += ' & !' + str(13 + i)
    fcouplecont += '\n'
    fcouplecont += 'q\n'
    fcouple.write(fcouplecont)
    fcouple.close()

    # generate input_nopbc.txt
    fpbcname = os.path.join(outpath, 'Share', 'input_nopbc.txt')
    fpbc = open(fpbcname, 'w')
    fpbc.write('1\n0\n')
    fpbc.close()

    # generate input_fit.txt
    ffitname = os.path.join(outpath, 'Share', 'input_fit.txt')
    ffit = open(ffitname, 'w')
    ffit.write('4\n0\n')
    ffit.close()



def main():
    parser = parseAgruments()
    args= vars(parser.parse_args())
    if len(args) < 6:
        parser.print_help()
    else:
        debug = args['debug']
    
    metadict = parseMetaFile(args['meta'], args['outpath'],debug)
    prodict, ligdict = getLigands(args['pdb'], debug)
    try:
        ff = metadict['ff']
        if debug:
            print("Looking for ions file at {}".format(os.path.join(args['data'], 'ions_'+ff+'.json')))
        ionsdict = json.load(open(os.path.join(args['data'], 'ions_'+ff+'.json')))
    except:
        print("There is not ions list found. Consider no type of ions")
        ionsdict = {}

    savedlig = separateLigands(prodict, ligdict, metadict, ionsdict, args['outpath'], args['data'], debug)
    processLigand(savedlig, args['outpath'], args['chimera'], debug)
    genInputGmx(args['outpath'], debug)

if __name__ == "__main__":
    main()
    