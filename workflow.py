import argparse as ap
from ast import arg
from curses import meta
from logging import root
import shutil
import sys
import os

def parseAgruments():
    parser = ap.ArgumentParser(description='Run MD simulation  \
         with several mutated proteins') 
    parser.add_argument('-w', '--workspace', type=str, nargs='?', \
        help='Workspace folder', \
            required=True)
    parser.add_argument('-i', '--input', type=str, nargs='?', \
        help='Input file contains parameters', \
        required=True)
    parser.add_argument('-f', '--fasta', type=str, nargs='?', \
        help='Fasta file of the wild type protein', \
        required=True)
    parser.add_argument('-m', '--mutation', type=str, nargs='?', \
        help='Mutation file contains list of muatation of interest', \
        required=True)
    parser.add_argument('-p', '--pdb', type=str, nargs='?', \
        help='Template model pdb file of the wildtype protein', \
        required=True)
    parser.add_argument('-a', '--alphafold', type=str, nargs='?', \
        help='Path to alphafold folder', \
        required=True)
    parser.add_argument('-c', '--docker', type=str, nargs='?', \
        help='Name of the docker image of workflow', \
        default='workflow')
    parser.add_argument('-l', '--acpype', type=str, nargs='?', \
        help='Name of the docker image of acpype', \
        default='acpype/acpype')
    parser.add_argument('-d', '--debug', help="Print log to debug or not", action='store_true')
    
    return parser

def genMutation(fasta, mutfile, outpath, debug=False):
    if debug:
        print('GENERATE MUTATION')
    command = ['python3', 'pythonScript/genMutation.py', '-f', fasta, '-m', mutfile, '-o', outpath]
    if debug:
        command.extend(['-d'])
    os.system(' '.join(command))
    if debug:
        print('\n\n')

def genShareFolder(metafile, pdbfile, rootpath, datafolder, dockername, debug=False):
    if debug:
        print('GENERATE SHARE FOLDER')
    # need to generate an sh file and then run a docker container to run it
    if not os.path.isdir(rootpath):
        os.mkdir(rootpath)
    # copy everything needed to run the script genShareFolder to rootpath
    # rootpath will be mounted to container under /tmp/workspace
    # copy data folder 
    datacont = os.path.join(rootpath, 'data')
    if os.path.isdir(datacont):
        if debug:
            print("Folder data already existed, delete the old one")
        shutil.rmtree(datacont)
    shutil.copytree(datafolder, datacont) 

    if debug:
        print("Finish copy data folder to {}".format(rootpath))

    # copy the script genShareFolder
    shutil.copy('pythonScript/genShareFolder.py', rootpath)
    if debug:
        print("Finish copy pythonScript/genShareFolder.py to {}".format(rootpath))

    # copy the pdbfile 
    shutil.copy(pdbfile, rootpath)
    if debug:
        print("Finish copy {} to {}".format(pdbfile, rootpath))

    # copy metafile 
    shutil.copy(metafile, rootpath)
    if debug:
        print("Finish copy {} to {}".format(metafile, rootpath))

    metaname = metafile.split('/')[-1]
    pdbname = pdbfile.split('/')[-1]


    # create run_genShare.py so the container can call it
    runfile = open(os.path.join(rootpath, 'run_genShare.sh'), 'w')
    runfile.write('#!/bin/bash\n')
    runfile.write('WS=/tmp/workspace/ \n')
    runfile.write('cd $WS\n')
    command = 'python3 genShareFolder.py -f ' + metaname + ' -p ' + \
        pdbname + ' -o /tmp/workspace/md ' + ' -i data -c /tmp/chimera/bin/chimera '
    if debug:
        command += ' -d\n'
    else:
        command += '\n'

    runfile.write(command)
    runfile.close()

    # now run the docker container to chmod the run_genShare.sh file
    # mount and bind rootpath to /tmp/workspace
    comchmod = 'docker run --gpus all --rm --user $(id -u):$(id -g) --mount src=' 
    comchmod += rootpath + ','
    comchmod += 'target=/tmp/workspace,type=bind '
    comchmod += dockername + ' '
    comchmod += 'chmod +x /tmp/workspace/run_genShare.sh'
    if debug:
        print(comchmod)
    os.system(comchmod)

    comrun = 'docker run --gpus all --rm --user $(id -u):$(id -g) --mount src=' 
    comrun += rootpath + ','
    comrun += 'target=/tmp/workspace,type=bind '
    comrun += dockername + ' '
    comrun += '/tmp/workspace/run_genShare.sh'
    if debug:
        print(comrun)
    os.system(comrun)

    if debug:
        print('\n\n')

    return os.path.join(rootpath, 'md', 'Share')



def parseMetaFile(metapath, debug=False):
    meta = {'type': 'homo', 'ff': 'amber', 'num_chain': 1, 'lig': '', \
        'ff_name': "amber99", 'water': 'tip3p', 'bt': "dodecahedron", \
        'p_name': "NA", 'n_name': "CL", 'sim_step': 2000000, \
        'data_dir': '/data/genetic_databases/'} # default meta data 
    if not os.path.isfile(metapath):
        print("The path of the input file {} is not correct".format(metapath))
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
    
    return meta

def predictBatch(metadict, fastapath, outpath, alphapath, datapath, debug=False):
    if debug:
        print('PREDICT STRUCTURE FOR MUTATED PROTEIN')

    model_preset = 'monomer'

    if metadict['type'] == 'homo' or metadict['type'] == 'complex':
        model_preset = 'multimer'


    if not os.path.isdir(outpath):
        if debug:
            print("Creating output directory to store AlphaFold result")
        os.mkdir(outpath)
    command = ['python3', 'pythonScript/predictBatch.py', '-f', fastapath, \
        '-o', outpath, '-a', alphapath, '-s', datapath, '-m', model_preset]
    if debug:
        command.extend (['-d'])
    
    os.system(' '.join(command))

    if debug:
        print('\n\n')

    
def runAcpype(sharepath, acpype, debug=False): 
    if debug:
        print('RUN ACPYPE WITH LIGANDS')
    if not os.path.isdir(sharepath):
        print("Cannot find Share directory")
        exit()
    if not os.path.isfile(os.path.join(sharepath, 'run_acpype.sh')):
        print("Cannot find run_acpype.sh in Share folder")
        exit()
    shareabs = os.path.abspath(sharepath)
    bind = 'src=' + shareabs + ',target=/tmp/workspace,type=bind'
    comchmod = ['docker', 'run', '--gpus', 'all', '--env', 'NVIDIA_DISABLE_REQUIRE=1', \
        '--rm', '--user', '$(id -u):$(id -g)', '--mount', \
        bind, acpype, 'chmod', '+x', '/tmp/workspace/run_acpype.sh']
    if debug:
        print(' '.join(comchmod))
    os.system(' '.join(comchmod))

    comrun =  ['docker', 'run', '--gpus', 'all', '--env', 'NVIDIA_DISABLE_REQUIRE=1', \
        '--rm', '--user', '$(id -u):$(id -g)', '--mount', \
        bind, acpype, '/tmp/workspace/run_acpype.sh']
    if debug:
        print(' '.join(comrun))
    os.system(' '.join(comrun))

    if debug:
        print('\n\n')

def processLigand(sharepath, debug=False):
    if debug:
        print('POST PROCESSING LIGAND AFTER ACPYPE')
        print("Current folder is {}".format(sharepath))

    subfolders = os.scandir(sharepath)

    for subfolder in subfolders:
        if os.path.isdir(subfolder) and '.acpype' in subfolder.name and subfolder.name[0] != '.':
            fullpath = os.path.abspath(subfolder)
            command = ['python3', 'pythonScript/processLigandsAcpype.py', '-f', \
                fullpath]
            if debug:
                command.extend(['-d'])
            os.system(' '.join(command))

    if debug:
        print('\n\n')

def alignStructure(rootpath, pdbfile, docker, debug=False):
    if debug:
        print('ALIGN PREDICTED STRUCTURE WITH TEMPLATE MODEL')
    # bind rootpath to docker container and then use chimera from the container to align structure

    # copy the script alignStructure.py to rootpath
    shutil.copy('pythonScript/alignStructure.py', rootpath)
    if debug:
        print("Finish copy pythonScript/alignStructure.py to {}".format(rootpath))
        
    # copy the pdbfile 
    shutil.copy(pdbfile, rootpath)
    if debug: 
        print("Finish copy {} to {}".format(pdbfile, rootpath))
    pdbname = pdbfile.split('/')[-1]

    # generate run_align.sh script to call alignStructure.py 
    script = '#!/bin/bash\n'
    script += 'WS=/tmp/workspace\n'
    script += 'cd $WS\n'
    script += 'python3 alignStructure.py -f alres -o md -c /tmp/chimera/bin/chimera -p ' + pdbname
    if debug:
        script += ' -d\n'
    else:
        script += '\n'

    runfile = open(os.path.join(rootpath, 'run_align.sh'), 'w')
    runfile.write(script)
    runfile.close()

    bind = 'src=' + rootpath + ',target=/tmp/workspace,type=bind'

    comchmod = ['docker', 'run', '--gpus', 'all', '--rm', '--user', '$(id -u):$(id -g)', \
        '--mount', bind, docker, 'chmod', '+x', '/tmp/workspace/run_align.sh']
    if debug:
        print(' '.join(comchmod))
    os.system(' '.join(comchmod))

    comrun = ['docker', 'run', '--gpus', 'all', '--rm', '--user', '$(id -u):$(id -g)', \
        '--mount', bind, docker, '/tmp/workspace/run_align.sh']
    if debug:
        print(' '.join(comrun))
    os.system(' '.join(comrun))

    if debug:
        print('\n\n')

# image is name of the docker image
def runMDBatch(rootpath, image, debug=False):
    if debug:
        print('RUN MD SIMULATION AND BINDING FREE ENERGY CALCULATION')
    # mount rootpath to docker container and run md simulation inside the md subfolder
    # check md subfolder and Share folder inside md first 
    if not os.path.isdir(os.path.join(rootpath, 'md')):
        print("There is no 'md' folder inside {}".format(rootpath))
        exit()
    if not os.path.isdir(os.path.join(rootpath, 'md', 'Share')):
        print("There is no 'Share' folder inside {}".format(os.path.join(rootpath, 'md')))
        exit()

    #copy MDBatch.py to rootpath to run 
    shutil.copy('pythonScript/MDBatch.py', rootpath)
    if debug:
        print('Copy MDBatch.py to {}'.format(rootpath))

    # rootpath will be mounted to /tmp/workspace/ inside the container 
    # create file run_md.sh for container to call
    
    content = '#!/bin/bash\n'
    content += 'WS=/tmp/workspace/ \n'
    content += 'cd $WS \n'
    content += 'source /tmp/amber20/amber.sh \n'
    content += 'export GMX=/usr/local/gromacs \n'
    content += 'export PATH="$GMX/bin:$PATH" \n'
    content += 'python3 MDBatch.py -f md -g gmx_mpi -a /tmp/amber20/miniconda/bin/gmx_MMPBSA -m 3'
    if debug:
        content += ' -d \n'
    else:
        content += '\n'

    runfile = open(os.path.join(rootpath, 'run_MD.sh'), 'w')
    runfile.write(content)
    runfile.close()

    # now mount rootpath to docker image
    rootabs = os.path.abspath(rootpath) 
    bind = 'src=' + rootabs + ',target=/tmp/workspace,type=bind'
    comchmod = ['docker', 'run', '--gpus', 'all', '--rm', '--user', '$(id -u):$(id -g)', \
        '=--mount', bind, image, 'chmod', '+x', '/tmp/workspace/run_MD.sh']
    if debug:
        print(' '.join(comchmod))
    os.system(' '.join(comchmod))

    comrun = ['docker', 'run', '--gpus', 'all', '--rm', '--user', '$(id -u):$(id -g)', \
        '--mount', bind, image, '/tmp/workspace/run_MD.sh']
    if debug:
        print(' '.join(comrun))
    os.system(' '.join(comrun))
    
    if debug:
        print('\n\n')

# delete files and folder generated during running time
def cleanup(rootpath, debug=False):
    if debug:
        print('CLEAN UP')
    if os.path.isfile(os.path.join(rootpath, 'alignStructure.py')):
        if debug:
            print("Delete file {}".format('alignStructure.py'))
        os.remove(os.path.join(rootpath, 'alignStructure.py'))
    
    if os.path.isfile(os.path.join(rootpath, 'genShareFolder.py')):
        if debug:
            print("Delete file {}".format('genShareFolder.py'))
        os.remove(os.path.join(rootpath, 'genShareFolder.py'))
    
    if os.path.isfile(os.path.join(rootpath, 'run_align.sh')):
        if debug:
            print("Delete file {}".format('run_align.sh'))
        os.remove(os.path.join(rootpath, 'run_align.sh'))

    if os.path.isfile(os.path.join(rootpath, 'run_MD.sh')):
        if debug:
            print("Delete file {}".format('run_MD.sh'))
        os.remove(os.path.join(rootpath, 'run_MD.sh'))
    
    if os.path.isfile(os.path.join(rootpath, 'MDBatch.py')):
        if debug:
            print("Delete file {}".format('MDBatch.py'))
        os.remove(os.path.join(rootpath, 'MDBatch.py'))

    if os.path.isfile(os.path.join(rootpath, 'run_genShare.sh')):
        if debug:
            print("Delete file {}".format('run_genShare.sh'))
        os.remove(os.path.join(rootpath, 'run_genShare.sh'))
    


    

def main():
    parser = parseAgruments()
    args= vars(parser.parse_args())
    if len(args) < 9:
        parser.print_help()
    else:
        debug = args['debug']
        jobname = args['pdb'].split('/')[-1].replace('.pdb', '')
        if debug:
            print('\n\n')
            print("WORKING WITH {}".format(jobname.upper()))
        scriptpath = os.getcwd() 
        rootpath = os.path.abspath(os.path.join(args['workspace'], jobname))

        if debug:
            print('Root directory is {}'.format(rootpath))
        if not os.path.isdir(rootpath):
            print("Creating {}".format(rootpath))
            os.mkdir(rootpath)

        fastapath = os.path.join(args['workspace'], jobname, 'fastas')
        if not os.path.isdir(fastapath):
            print("Creating {} to save fasta files of mutations".format(fastapath))
            os.mkdir(fastapath)
        
        alphares = os.path.join(rootpath, 'alres')
        if not os.path.isdir(alphares):
            print("Creating {} to save AlphaFold2's results".format(alphares))
            os.mkdir(alphares)

        metadic = parseMetaFile(args['input'], debug)

        if debug:
            print('\n\n')

        genMutation(args['fasta'], args['mutation'], fastapath, debug)
        predictBatch(metadic, fastapath, alphares, args['alphafold'], metadic['data_dir'], debug)
        sharepath = genShareFolder(args['input'], args['pdb'], rootpath, \
            'pythonScript/data', args['docker'], debug)

        sharepath = os.path.join(rootpath, 'md', 'Share')
        runAcpype(sharepath, args['acpype'], debug)
        processLigand(sharepath, debug)
        alignStructure(rootpath, args['pdb'], args['docker'], debug)
        runMDBatch(rootpath, args['docker'], debug)
        cleanup(rootpath, debug)
        
if __name__ == "__main__":
    main()
