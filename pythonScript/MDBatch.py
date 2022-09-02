import argparse as ap
from logging import exception
import sys
import os
import fileinput


def parseAgruments():
    parser = ap.ArgumentParser(description='Run MD simulation  \
         with several mutated proteins') 
    parser.add_argument('-f', '--folder', type=str, nargs='?', \
        help='Path of the folder contains a "Share" folder and also subfolders of mutation proteins, \
            inside each subfolder there is a pdb file of protein, \
            this subfolder is also the folder to run md for the protein inside \
            all the output will be generated here', \
            required=True)
    parser.add_argument('-g', '--gmx', type=str, nargs='?', default='gmx', \
         help='Command to call gmx') 
    parser.add_argument('-a', '--mmpbsa', type=str, nargs='?', default='gmx_mmpbsa', \
         help='Command to call gmx_mmpbsa') 
    parser.add_argument('-m', '--mode', type=int, nargs='?', default=3, \
        help='Mode: 1 - compute MD simulation only, \
            2 - compute MMPBSA only, \
            3 - compute both MD and MMPBSA')
    parser.add_argument('-d', '--debug', help="Print log to debug or not", action='store_true')
    
    return parser

def readLigandFile(path, debug=False):
    lfile = open(path, 'r')
    lines = lfile.readlines() 
    try: 
        nspecies = int(lines[0])
    except:
        print("Number of species is not an integer, receive from input file: {}".format(lines[0]))
        exit()
    # keep track the index of the next line to be read, 
    # omit the seperate line ====== at line_ind = 1
    line_ind = 2
    list_species = []
    if debug: 
        print("Number of ligand species: {}".format(nspecies))
    for i in range(nspecies):
        species = {}
        species['name'] = lines[line_ind].strip()
        line_ind += 1 # keep track of the next line to be read
        if debug:
            print("Work with ligand type {}, named {}".format(i, species['name']))

        try:
            nligands = int(lines[line_ind])
        except:
            print("Number of ligands of type {} is not an integer, receive from input file: {}".format(i, lines[line_ind]))
            exit()

        line_ind += 1 # keep track of the index of next line to be read
        species['nligands'] = nligands
        if debug:
            print("There are {} ligands of type {}".format(nligands, species['name']))
        list_gro = []
        for j in range(nligands):
            list_gro.append (lines[line_ind].strip())
            if '.gro' not in lines[line_ind]:
                print("The format of the coordinate file is not .gro, get: {}".format(lines[line_ind]))
            line_ind += 1 # keep track of the next line to be read 
        species['gros'] = list_gro
        species['itp'] = lines[line_ind].strip()
        if '.itp' not in species['itp']:
            print('The format of ligand topology file is not .itp, get: {}'.format(species['itp']))
            exit() 
        line_ind += 1 #keep track of the next line to be read 

        species['prm'] = lines[line_ind].strip()
        if '.prm' not in species['prm']:
            print('The format of ligand parameter file is not .itp, get: {}'.format(species['itp']))
            exit() 
        line_ind += 1 #keep track of the next line to be read 

        species['posre'] = lines[line_ind].strip()
        if '.itp' not in species['posre']:
            print('The format of ligand position restraint file is not .itp, get: {}'.format(species['itp']))
            exit() 
        line_ind += 2 # omit the separate line ====== between species
        list_species.append(species)
    if debug:
        print(list_species) 
    return list_species

def readConfigFile(path, debug=False):
    cfile = open(path, 'r')
    conf = {'pdb2gmx': '', \
            'editconf':'', \
            'genion': '' , \
            'grompp': ''}
    lines = cfile.readlines()
    for line in lines:
        if 'pdb2gmx' in line.strip().lower():
            conf['pdb2gmx'] = line.strip().split(':')[1]
        if 'editconf' in line.strip().lower():
            conf['editconf'] = line.strip().split(':')[1]
        if 'genion' in line.strip().lower():
            conf['genion'] = line.strip().split(':')[1]
        if 'grompp' in line.strip().lower():
            conf['grompp'] = line.strip().split(':')[1]
    if debug:
        print(conf)
    return conf

def modifyTopol(path, list_species, debug=False):
    for line in fileinput.input(path,inplace=1):
        if 'forcefield.itp' in line:
            new_lines = '\n; Include ligand parameters .prm files \n'
            for species in list_species:
                new_lines += '#include "'
                new_lines += species['prm'] + '"\n'
            line = line.replace(line, line+new_lines)
        
        if 'Include water topology' in line:
            new_lines = '; Include ligand topology and position restraint\n'
            for species in list_species:
                new_lines += '#include "'
                new_lines += species['itp'] + '"\n'
                new_lines += "#ifdef POSRES\n"
                new_lines += '#include "'
                new_lines += species['posre'] + '"\n'
                new_lines +=  "#endif\n\n"

            new_lines += '; Include water topology\n'
            line = line.replace(line, new_lines)

        # if '#mols' in line:
        #     new_lines = ''
        #     for species in list_species:
        #         new_lines += species['name']
        #         new_lines += '\t\t\t\t\t'
        #         new_lines += str (species['nligands'])
        #         new_lines += '\n'
        #     line = line.replace(line, line+new_lines)
        sys.stdout.write(line)

    file = open(path, 'a')
    for species in list_species:
        new_line = species['name'] + '\t\t\t\t\t' + str(species['nligands']) + '\n'
        file.write(new_line)

def modifyGro(path, list_species, debug=False):
    fpro = open(path, 'r')
    prolines = fpro.readlines()
    numAtomPro = int(prolines[1])
    if debug:
        print("Number of atoms of the protein: {}".format(numAtomPro))
    line_to_add = []
    for species in list_species:
        list_gros = species['gros']
        for gro in list_gros:
            # read file cooridnate of each ligand 
            flig = open(gro, 'r')
            liglines = flig.readlines()
            numAtomLig = int(liglines[1])
            if debug:
                print("Number of atoms of the ligands: {}".format(numAtomLig))
            line_to_add.extend(liglines[2:2+numAtomLig])
    if debug:
        print("Adding {} new lines".format(len(line_to_add)))
        # print(line_to_add)
    comlines = [prolines[0]]
    comlines.extend([str(numAtomPro + len(line_to_add)) + '\n']) # new number of atoms 
    comlines.extend(prolines[2:2 + numAtomPro])
    comlines.extend(line_to_add)
    comlines.extend([prolines[len(prolines)-1]])
    # print(comlines)
    with open('com.gro', 'w') as fw:
        fw.write(''.join(comlines))
    
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
        'aid': id, 'atom': atom, 'name': name, \
        'chain': chain, 'rid': resid, \
        'line': line}
    return linedict

def runPdb2gmx(rootpath, gmx, pdbfile, config, debug=False):
    # add ions here if needed 
    if os.path.isfile(os.path.join(rootpath, "Share", 'ions.txt')):
        if debug:
            print("There are some ions to add to the .pdb file")
        fions = open(os.path.join(rootpath, 'Share', 'ions.txt'), 'r')
        ionslines = fions.readlines()
        fpdb = open(pdbfile, 'r')
        pdblines = fpdb.readlines()
        fpdb.close()

        # find the position to insert 
        pos = len(pdblines) - 1
        ister = False
        for i in range(len(pdblines) - 1, -1, -1):
            if 'TER' in pdblines[i]:
                ister = True
            if 'ATOM' in pdblines[i] or 'TER' in pdblines[i]:
                pos = i
                if pdblines[i][-1] != '\n':
                    pdblines[i] += '\n'
                break 

        if not ister: # in case if there is no TER line, add TER line
            lastatom = pdblines[pos]
            linedict = parseLine(lastatom)
            terline = "TER  "
            for j in range(6 - len(str(linedict['aid']+1))):
                terline += ' '
            terline += str(linedict['aid']+1)
            terline += '      '
            terline += linedict['name']
            terline += ' '
            terline += linedict['chain']
            for j in range(4 - len(str(linedict['rid']))):
                terline += ' '
            terline += str(linedict['rid'])
            terline += '\n'
            pdblines.insert(pos+1, terline)
            pos += 1
            
        lastline = pdblines[pos]
        lastlinedict = parseLine(lastline)
        lastid = lastlinedict['aid']

        addedline = 0
        for ionline in ionslines:
            if 'HETATM' in ionline: 
                # need to modify the id of the ions to be consistent with atom lines
                iondict = parseLine(ionline)
                ionid = iondict['aid']
                if len(str(ionid)) == len(str(lastid + 1)):
                    ionline = ionline.replace(str(ionid), str(lastid+1), 1)
                elif len(str(ionid)) < len(str(lastid+1)):
                    todelete = ''
                    for j in range(len(str(lastid+1)) - len(str(ionid))):
                        todelete += ' '
                    todelete += str(ionid)
                    ionline = ionline.replace(todelete, str(lastid+1), 1)
                else: # len(str(ionid)) > len(str(lastid+1))
                    toadd = ''
                    for j in range(len(str(ionid))-len(str(lastid+1))):
                        toadd += ' '
                    toadd += str(lastid+1)
                    ionline = ionline.replace(str(ionid), toadd, 1)
                pdblines.insert(pos + addedline + 1, ionline)
                addedline += 1
                lastid += 1

        fpdb = open(pdbfile, 'w')
        fpdb.write(''.join(pdblines))
        fpdb.close() 


    optPdb2gmx = config['pdb2gmx'].split()
    command = [gmx, 'pdb2gmx', '-f', pdbfile, '-o', 'protein.gro'] 
    command.extend(optPdb2gmx)
    if debug:
        print("Command to call:")
        print(" ".join(command))
    os.system(' '.join(command)) # command this for test 
    
def runSolvation(gmx, config, debug=False):
    if config:
        opts = config['editconf'].split()
        if debug:
            print("Solvate the system with options: {}".format(config))
        edifconf = [gmx, 'editconf', '-f', 'com.gro', '-o', 'box.gro']
        edifconf.extend(opts)
    else:
        edifconf = [gmx, 'editconf', '-f', 'com.gro', \
            '-o', 'box.gro', '-bt', 'dodecahedron', '-d', '1.0']
    if debug:
        print("Command to editconf the system")
    os.system(" ".join(edifconf))

    solv = [gmx, 'solvate', '-cp', 'box.gro', '-cs', 'spc216.gro',  \
        '-p', 'topol.top', '-o',  'solv.gro']
    os.system(" ".join(solv))

def runNeutralize(gmx, mdp, config, input, debug=False):
    if not os.path.isfile(mdp):
        print("There is no mdp file for adding ions")
        exit()

    if config['grompp']:
        grompp = [gmx, 'grompp', '-f', mdp, '-c', 'solv.gro', \
            '-p', 'topol.top', '-o', 'ions.tpr']
        grompp.extend(config['grompp'].split())
    else:
        grompp = [gmx, 'grompp', '-f', mdp, '-c', 'solv.gro', \
            '-p', 'topol.top', '-o', 'ions.tpr']
    os.system(" ".join(grompp))

    if config['genion']:
        genion = [gmx, 'genion',  '-s', 'ions.tpr', '-o', 'solv_ions.gro', \
            '-np ', '1','-p', 'topol.top', '-neutral']
        genion.extend(config['genion'].split())
    else:
        genion = [gmx, 'genion',  '-s', 'ions.tpr', '-o', 'solv_ions.gro', \
            '-np ', '1', '-p', 'topol.top', '-pname', 'NA', '-nname', 'CL' '-neutral']
    
    if not os.path.isfile(input):
        print("There is no input file for genion, input by stdin")
    else:
        genion.extend(['<', input])

    os.system(" ".join(genion))


def runEM(gmx, mdp, config, debug=False):
    if not os.path.isfile(mdp):
        print("There is not mdp file for energy minimization")
        exit()
    gromm = [gmx, 'grompp', '-f', mdp, '-c', \
        'solv_ions.gro', '-p', 'topol.top', '-o', 'em.tpr']
    if config['grompp']:
        gromm.extend(config['grompp'].split())
    os.system(' '.join(gromm))

    em = [gmx, 'mdrun', '-v', '-deffnm', 'em']
    os.system(' '.join(em))


def runCoupling(gmx, input, debug = False):
    couple = [gmx, 'make_ndx', '-f', 'em.gro', '-o', 'index.ndx']

    if not os.path.isfile(input):
        print("There is no input file for coupling, input by stdin")
    else:
        couple.extend(['<', input])

    os.system(' '.join(couple))

def runEquilMD(gmx, nvtmdp, nptmdp, mdmdp, config, debug=False):
    if not os.path.isfile(nvtmdp):
        print("No mdp file for nvt equilibrium")
        exit()
    if not os.path.isfile(nptmdp):
        print("No mdp file for npt equilibrium")
        exit()
    if not os.path.isfile(mdmdp):
        print("No mdp file for md")
        exit()

    nvtpp = [gmx, 'grompp', '-f', nvtmdp, '-c', 'em.gro',\
        '-r', 'em.gro', '-p', 'topol.top', '-n', 'index.ndx', '-o', 'nvt.tpr']
    if config['grompp']:
        nvtpp.extend(config['grompp'].split())
    os.system(' '.join(nvtpp))

    nvtmd = [gmx, 'mdrun', '-deffnm', 'nvt']
    os.system(" ".join(nvtmd))

    nptpp = [gmx, 'grompp', '-f', nptmdp, '-c', 'nvt.gro', \
        '-t', 'nvt.cpt', '-r', 'nvt.gro', '-p', 'topol.top', '-n', 'index.ndx', '-o', 'npt.tpr']
    if config['grompp']:
        nptpp.extend(config['grompp'].split())
    os.system(' '.join(nptpp))

    nptmd = [gmx, 'mdrun', '-deffnm', 'npt']
    os.system(" ".join(nptmd))

    mdpp = [gmx, 'grompp', '-f', mdmdp, '-c', 'npt.gro', '-t', \
        'npt.cpt', '-p', 'topol.top', '-n', 'index.ndx', '-o', 'md.tpr']
    if config['grompp']:
        mdpp.extend(config['grompp'].split())
    os.system(" ".join(mdpp))

    md = [gmx, 'mdrun', '-deffnm', 'md']
    os.system(" ".join(md))

def runNopbc(gmx, input, debug=False):
    nopbc = [gmx, 'trjconv', '-s', 'md.tpr', '-f', 'md.xtc', '-o', \
        'md_nopbc.xtc', '-center', '-pbc', 'mol', '-ur', 'compact']

    if not os.path.isfile(input):
        print("There is no input file for coupling, input by stdin")
    else:
        nopbc.extend(['<', input])    

    os.system(' '.join(nopbc))

def runFit(gmx, input, debug=False):
    fit = [gmx, 'trjconv', '-s', 'md.tpr', '-f', 'md_nopbc.xtc', \
        '-o', 'md_fit.xtc', '-fit', 'rot+trans']

    if not os.path.isfile(input):
        print("There is no input file for coupling, input by stdin")
    else:
        fit.extend(['<', input])    
    os.system(' '.join(fit))

def processMD(rootpath, gmx, list_species, conf, debug=False):
    # Reading config files 
    ionmdp = os.path.join(rootpath, 'Share', 'ions.mdp')
    emmdp = os.path.join(rootpath, "Share", "em.mdp")
    nvtmdp = os.path.join(rootpath, "Share", "nvt.mdp")
    nptmdp = os.path.join(rootpath, "Share", "npt.mdp")
    mdmdp = os.path.join(rootpath, "Share", 'md.mdp') 

    genioninput = os.path.join(rootpath, "Share", 'input_genion.txt')
    couplinginput = os.path.join(rootpath, "Share", 'input_coupling.txt')
    nopbcinput = os.path.join(rootpath, "Share", 'input_nopbc.txt')
    fitinput = os.path.join(rootpath, "Share", 'input_fit.txt')


    subfolders = os.scandir(rootpath)
    
    for subfolder in subfolders:
        if subfolder.name[0] != '.' and subfolder.name != 'Share' and os.path.isdir(subfolder):
            if debug:
                print("Run MD in folder {}".format(subfolder.name))
            files = os.scandir(subfolder)
            for file in files:
                if '.pdb' in file.name:
                    if debug:
                        print('Working with file {}'.format(file.name))
                    # change to current directory
                    os.chdir(subfolder)
                    if debug:
                        print("Current directory is: {}".format(os.getcwd()))

                    # Run pdb2gmx
                    runPdb2gmx(rootpath, gmx, file.name, conf, debug)
                    
                    # add ligand into the com.gro files 
                    modifyTopol('topol.top', list_species, debug) # command this for test 
                    modifyGro('protein.gro', list_species, debug)
                    # Solvate the system 
                    runSolvation(gmx, conf, debug)
                    runNeutralize(gmx, ionmdp, conf, genioninput, debug)
                    runEM(gmx,emmdp, conf, debug)
                    runCoupling(gmx, couplinginput, debug)
                    runEquilMD(gmx, nvtmdp, nptmdp, mdmdp, conf, debug)
                    runNopbc(gmx, nopbcinput, debug)
                    runFit(gmx, fitinput, debug)
                    

def processMMPBSA(rootpath, mmpbsa, debug=False):
    input = os.path.join(rootpath, "Share", 'mmpbsa.in')
    if not os.path.isfile(input):
        print("There is no input file to run MMPBSA")
        exit()

    subfolders = os.scandir(rootpath)
    
    for subfolder in subfolders:
        if subfolder.name[0] != '.' and subfolder.name != 'Share' and os.path.isdir(subfolder):
            os.chdir(subfolder)
            if debug:
                print("Run MMPBSA in folder {}".format(subfolder.name))
                os.system('pwd')
            
            mm = [mmpbsa, 'MPI', '-O', '-i', input, '-cs', 'md.tpr', '-ci', 'index.ndx', \
                '-cg', '1', '13', '-ct', 'md_fit.xtc', '-cp', 'topol.top', '-nogui']
            os.system(" ".join(mm))

def main():
    parser = parseAgruments()
    args= vars(parser.parse_args())
    if len(args) < 5:
        parser.print_help()
    else:
        debug = args['debug']
        if debug:
            print('========================Start processing========================')
            print('Parsed enough arguments, processing')
        
        # get absolute path of all the input directory 
        absin = os.path.abspath(args['folder'])
        if debug:
            print("The absolute paths:")
            print(absin)
        
        # check if the paths are correct 
        if not os.path.isdir(absin):
            print("The input folder is not correct")
            exit()
        
        # check existance of the ligands.txt file
        if not os.path.isfile(os.path.join(absin, 'Share', "ligands.txt")):
            print("There is no ligand file at {}, check the path".format(os.path.join(absin, 'Share', "ligands.txt")))
            exit()
        else: 
            if debug:
                print("The ligand file is at {}".format(os.path.join(absin, 'Share')))
        
        # check existance of the config_pdb2gmx.txt file
        if not os.path.isfile(os.path.join(absin, 'Share', "config.txt")):
            print("There is no config file at {}, check the path".format(os.path.join(absin, 'Share', "config.txt")))
            exit()
        else: 
            if debug:
                print("The config file is at {}".format(os.path.join(absin, 'Share')))


    # read ligands.txt file 
    pathlfile = os.path.join(absin, 'Share', 'ligands.txt')
    pathcfile = os.path.join(absin, 'Share', "config.txt")
    list_species = readLigandFile(pathlfile, debug=debug)
    conf = readConfigFile(pathcfile, debug)
    if not list_species:
        print("Empty ligand information")
        exit()
    if debug:
        if args['mode'] == 1:
            print("Compute MD simulation only")
        if args['mode'] == 2:
            print("Compute MMPBSA only")
        if args['mode'] == 3:
            print("Compute both MD simulation and MMPBSA")

    if args['mode'] == 1 or args['mode'] == 3:
        # print("Run MD")
        processMD(absin, args['gmx'], list_species, conf, debug)
    if args['mode'] == 2 or args['mode'] == 3:
        processMMPBSA(absin, args['mmpbsa'], debug)
        # print("Run MMPBA")

                
if __name__ == "__main__":
    main()
    