import argparse as ap
from genericpath import isdir
import os
from genShareFolder import parseLine
from genShareFolder import getLigands
from genShareFolder import separateLigands
import json
from MDBatch import processMD
from MDBatch import processMMPBSA
import shutil
from genShareFolder import genInputGmx

def parseAgruments():
    parser = ap.ArgumentParser(description='Run MD simulation  \
         with several mutated proteins') 
    parser.add_argument('-w', '--workspace', type=str, nargs='?', \
        help='Absolute path to the workspace folder which will be mounted to the docker container', \
            required=True)
    parser.add_argument('-o', '--outfolder', type=str, nargs='?', \
        help='Output folder, path starts from workspace', \
            required=True)
    parser.add_argument('-i', '--input', type=str, nargs='?', \
        help='Input file contains parameters, path starts from workspace', \
        required=True)
    parser.add_argument('-p', '--pdbFolder', type=str, nargs='?', \
        help='Folder contains pdb files of the experimental protein to align with, path start from the workspace', \
        required=True)
    parser.add_argument('-a', '--alphaPredicted', type=str, nargs='?', \
        help='Pdb file predicted by alphafold', required=True)
    parser.add_argument('-f', '--data', type=str, nargs='?', \
        help='Path to data folder, from the workspace', required=True)
    parser.add_argument('-c', '--docker', type=str, nargs='?', \
        help='Name of the docker image of workflow', \
        default='workflow')
    parser.add_argument('-l', '--acpype', type=str, nargs='?', \
        help='Name of the docker image of acpype', \
        default='acpype/acpype')
    parser.add_argument('-d', '--debug', help="Print log to debug or not", action='store_true')
    
    return parser

def parsePdbFile(pdbfile, debug=False):
    return getLigands(pdbfile, debug)


def addHydrogen(workspace, sharefolder, docker, debug=False): # path of share folder is started from workspace folder 
    print("=========START ADDING HYDROGEN===========")
    acpypefile = '#!/bin/sh\n'
    acpypefile += 'WS=/tmp/workspace/\n' # starting folder of docker container
    acpypefile += 'cd $WS\n'

    lfiles = os.listdir(os.path.join(workspace, sharefolder))
    for file in lfiles:
        if '.pdb' in file:
            if debug:
                print("Adding hydrogen to the ligand {}".format(file))
            ligpath = os.path.join('/tmp/workspace', sharefolder, file)
            ligruncont = '/tmp/miniconda3/bin/obabel ' + ligpath + ' -O ' + ligpath + ' -h '
            ligrunfile = open(os.path.join(workspace, sharefolder, 'runAddHydrogen.sh'),'w')
            ligrunfile.write(ligruncont)
            ligrunfile.close()

            comrun = ['docker', 'run', '--env', 'NVIDIA_DISABLE_REQUIRE=1', '--rm', '--user',\
                '$(id -u):$(id -g)', '--mount', 'src=/home/t.pham1/workspace,target=/tmp/workspace,type=bind',  \
                docker, os.path.join('/tmp/workspace', sharefolder, 'runAddHydrogen.sh')] 
            if debug:
                print(' '.join(comrun))

            os.system(' '.join(comrun))
            
            os.remove(os.path.join(workspace, sharefolder, 'runAddHydrogen.sh'))

            acpypefile += 'acpype -f -i ' + file + ' -o gmx -d\n'

    facname = os.path.join(workspace, sharefolder, "run_acpype.sh")
    fac = open(facname, 'w')
    fac.write(acpypefile)
    fac.close()
    print("=========DONE ADD HYDROGEN===========")



def genProteinFolder(workspace, pdbfolder, outfolder, datafolder, meta, ions, acpype, docker, debug=False):
    if not os.path.isdir(os.path.join(workspace, pdbfolder)):
        print("Cannot find pbd folder")
        exit()
    if not os.path.isdir(os.path.join(workspace, outfolder)):
        os.mkdir(os.path.join(workspace, outfolder))
    lfiles = os.listdir(os.path.join(workspace, pdbfolder))

    for file in lfiles:
        if '.pdb' in file:
            if debug:
                print(file)
            protein = file.replace('.pdb', '')
            fullfile = os.path.join(workspace, pdbfolder, file)
            pros, ligs = parsePdbFile(fullfile, debug)
            if debug:
                print(pros.keys())
                print(ligs.keys())
            # create a folder for output 
            if not os.path.isdir(os.path.join(workspace, outfolder, protein)):
                os.mkdir(os.path.join(workspace, outfolder, protein))
            if not os.path.isdir(os.path.join(workspace, outfolder, protein, protein)):
                os.mkdir(os.path.join(workspace, outfolder, protein, protein))

            if meta['type'] == 'mono':
                key = list(pros.keys())[0] # take only the first one 
                
                proname = os.path.join(workspace, outfolder, protein, protein, 'protein.pdb')
                profile = open(proname, 'w')
                for line in pros[key]:
                    profile.write(line)
                profile.close()

                separateLigands(pros, ligs, meta, ions, os.path.join(workspace, outfolder, protein), \
                    datafolder, debug)
            
            addHydrogen(workspace, os.path.join(outfolder, protein, 'Share'), docker, debug)
            runAcpype(workspace, os.path.join(outfolder, protein, 'Share'), acpype, debug)
            processLigand(workspace, os.path.join(workspace, outfolder, protein, 'Share'), debug) 
            print("=========GENERATE INPUT FOR GMX=========")
            genInputGmx(os.path.join(workspace, outfolder, protein), debug)


def parseMetaFile(workspace, outfolder, inputfile, debug=False):
    meta = {'type': 'homo', 'ff': 'amber', 'num_chain': 1, 'lig': '', \
        'ff_name': "amber99", 'water': 'tip3p', 'bt': "dodecahedron", \
        'p_name': "NA", 'n_name': "CL", 'sim_step': 2000000, \
        'data_dir': '/data/genetic_databases/'} # default meta data 
    absinput = os.path.join(workspace, inputfile)
    if not os.path.isfile(absinput):
        print("The path of the input file {} is not correct".format(os.path.join(workspace, inputfile)))
        exit()
    file = open(os.path.join(workspace, inputfile), 'r')
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

    return meta

def runAcpype(workspace, sharepath, acpype, debug=False): 
    if debug:
        print('RUN ACPYPE WITH LIGANDS')
    if not os.path.isdir(os.path.join(workspace,sharepath)):
        print("Cannot find Share directory {}".format(sharepath))
        exit()
    if not os.path.isfile(os.path.join(workspace, sharepath, 'run_acpype.sh')):
        print("Cannot find run_acpype.sh in Share folder")
        exit()
    shareabs = os.path.abspath(os.path.join(workspace, sharepath))
    bind = 'src=' + shareabs + ',target=/tmp/workspace,type=bind'
    comchmod = ['docker', 'run', '--env', 'NVIDIA_DISABLE_REQUIRE=1', '--env', 'NVIDIA_DISABLE_REQUIRE=1', \
        '--rm', '--user', '$(id -u):$(id -g)', '--mount', \
        bind, acpype, 'chmod', '+x', '/tmp/workspace/run_acpype.sh']
    if debug:
        print(' '.join(comchmod))
    os.system(' '.join(comchmod))

    comrun =  ['docker', 'run', '--env', 'NVIDIA_DISABLE_REQUIRE=1', '--env', 'NVIDIA_DISABLE_REQUIRE=1', \
        '--rm', '--user', '$(id -u):$(id -g)', '--mount', \
        bind, acpype, '/tmp/workspace/run_acpype.sh']
    if debug:
        print(' '.join(comrun))
    os.system(' '.join(comrun))

    if debug:
        print('\n\n')

def processLigand(workspace, sharepath, debug=False):
    if debug:
        print('POST PROCESSING LIGAND AFTER ACPYPE')
        print("Current folder is {}".format(os.path.join(workspace,sharepath)))

    subfolders = os.scandir(os.path.join(workspace,sharepath))

    for subfolder in subfolders:
        if os.path.isdir(subfolder) and '.acpype' in subfolder.name and subfolder.name[0] != '.':
            fullpath = os.path.abspath(os.path.join(workspace, subfolder))
            command = ['python3', 'pythonScript/processLigandsAcpype.py', '-f', \
                fullpath]
            if debug:
                command.extend(['-d'])
            os.system(' '.join(command))

    if debug:
        print('\n\n')

def runMDBatch(workspace, outfolder, profolder, image, meta, debug=False):
    if debug:
        print('RUN MD SIMULATION AND BINDING FREE ENERGY CALCULATION')
    # mount rootpath to docker container and run md simulation inside the md subfolder
    # check md subfolder and Share folder inside md first 
    if not os.path.isdir(os.path.join(workspace, outfolder, profolder, 'Share')):
        print("There is no 'Share' folder inside \
            {}".format(os.path.join(workspace, outfolder, profolder, 'Share')))
        exit()

    #copy MDBatch.py to rootpath to run 
    shutil.copy('pythonScript/MDBatch.py', os.path.join(workspace, outfolder, profolder))
    if debug:
        print('Copy MDBatch.py to {}'.format(os.path.join(workspace, outfolder, profolder)))

    # profolder will be mounted to /tmp/workspace/ inside the container 
    print(profolder)
    # create file run_md.sh for container to call
    
    content = '#!/bin/bash\n'
    content += 'WS=/tmp/workspace/ \n'
    content += 'cd $WS \n'
    content += 'source /tmp/amber20/amber.sh \n'
    content += 'export GMX=/usr/local/gromacs \n'
    content += 'export PATH="$GMX/bin:$PATH" \n'
    content += 'python3 '+ profolder +'/MDBatch.py -f ' + profolder + \
        ' -g gmx_mpi -a /tmp/amber20/miniconda/bin/gmx_MMPBSA -m 3'
    if debug:
        content += ' -d \n'
    else:
        content += '\n'

    print(content)

    runfile = open(os.path.join(workspace, outfolder, profolder, 'run_MD.sh'), 'w')
    runfile.write(content)
    runfile.close()

    # write config file
    # generate config file from meta file 
    config = 'pdb2gmx: -ff ' + meta['ff_name'] + ' -water ' + meta['water'] + ' -ignh\n'
    config += 'edifconf: -bt ' + meta['bt']+ ' -d 1.0\n'
    config += 'genion: -pname ' + meta['p_name'] + ' -nname ' + meta['n_name'] + '\n'
    config += 'grompp: -maxwarn 50\n'

    if not os.path.isdir(os.path.join(workspace, outfolder, profolder, 'Share')):
        if debug:
            print("Creating 'Share' directory")
        os.mkdir(os.path.join(workspace, outfolder, profolder, 'Share'))
    fname = os.path.join(workspace, outfolder, profolder, 'Share','config.txt')
    fconfig = open(fname, 'w')
    fconfig.write(config)
    fconfig.close()

    # now mount rootpath to docker image
    rootabs = os.path.abspath(os.path.join(workspace, outfolder)) 
    bind = 'src=' + rootabs + ',target=/tmp/workspace,type=bind'
    comchmod = ['docker', 'run', '--env', 'NVIDIA_DISABLE_REQUIRE=1', '--rm', '--user', '$(id -u):$(id -g)', \
        '--mount', bind, image, 'chmod', '+x', '/tmp/workspace/' + profolder +'/run_MD.sh']
    if debug:
        print(' '.join(comchmod))
    os.system(' '.join(comchmod))

    comrun = ['docker', 'run', '--env', 'NVIDIA_DISABLE_REQUIRE=1', '--rm', '--user', '$(id -u):$(id -g)', \
        '--mount', bind, image, '/tmp/workspace/' + profolder +'/run_MD.sh']
    if debug:
        print(' '.join(comrun))
    os.system(' '.join(comrun))
    
    if debug:
        print('\n\n')


def main():
    parser = parseAgruments()
    args= vars(parser.parse_args())
    if len(args) < 8:
        parser.print_help()
    else:
        debug = args['debug']
        meta = parseMetaFile(args['workspace'], args['outfolder'], args['input'], debug)
        try:
            ff = meta['ff']
            if debug:
                print("Looking for ions file at {}".format(\
                    os.path.join(args['workspace'], args['data'], 'ions_'+ff+'.json')))
            ionsdict = json.load(open(os.path.join(\
                args['workspace'], args['data'], 'ions_'+ff+'.json')))
        except:
            print("There is not ions list found. Consider no type of ions")
            ionsdict = {}

        genProteinFolder(args['workspace'], args['pdbFolder'], args['outfolder'], \
            args['data'], meta, ionsdict, \
            args['acpype'], args['docker'], debug)
        
        fols = os.scandir(os.path.join(args['workspace'],args['outfolder']))
        for fol in fols:
            if os.path.isdir(os.path.join(args['workspace'], args['outfolder'], fol, 'Share')):
                # print(os.path.join(args['workspace'], args['outfolder'], fol, 'Share'))
                print("========START WITH {}==========".format(fol.name))
                runMDBatch(args['workspace'], args['outfolder'], fol.name, \
                    args['docker'], meta, debug)
                print("========DONE WITH {}==========".format(fol.name))

        


if __name__ == "__main__":
    main()
