import argparse as ap
from ast import parse
from datetime import date
from logging import exception
import numpy as np 
import os
from datetime import date
import subprocess


def parseAgruments():
    parser = ap.ArgumentParser(description='Predict batch of fasta files') 
    parser.add_argument('-f', '--folder', type=str, nargs='?', \
        help='Path of the folder contains fasta files', required=True)
    parser.add_argument('-o', '--outpath', type=str, nargs='?', \
        help='Path of output directory', required=True)
    parser.add_argument('-a', '--alphafold', type=str, nargs='?', \
         help='Path of alphafold directory (root directory)', required=True)

    today = date.today()
    todaystr = today.strftime("%Y-%m-%d")

    parser.add_argument('-t', '--templatedate', type=str, nargs='?', \
        help='Max template date for predicting protein structure, \
            default = ' + todaystr, default='2003-05-20')

    
    # parser.add_argument('-t', '--templatedate', type=str, nargs='?', \
    #     help='Max template date for predicting protein structure, \
    #         default = ' + todaystr, default=todaystr)
    parser.add_argument('-s', '--datadir', type=str, nargs='?', \
        help='Directory to the data folder, default = /data/genetic_databases/', default='/data/genetic_databases/')
    parser.add_argument('-m', '--modelpreset', type=str, nargs='?', \
        help="Model preset for alphafold, default = monomer", default='monomer')
    parser.add_argument('-r', '--relax', type=str, nargs = '?', \
        help='Run relaxation or not, default is true', default='true')
    parser.add_argument('-d', '--debug', help="Print log to debug or not", action='store_true')
    return parser
            

def main():
    parser = parseAgruments()
    args= vars(parser.parse_args())
    if len(args) < 2:
        parser.print_help()
    else:
        debug = args['debug']
        if debug:
            print('Parsed enough arguments, processing')
        
        # get absolute path of all the input directory 
        absin = os.path.abspath(args['folder'])
        absal = os.path.abspath(args['alphafold'])
        absout = os.path.abspath(args['outpath']) 
        if debug:
            print("The absolute paths:")
            print(absin)
            print(absal)
            print(absout)
        
        # check if the paths are correct 
        if not os.path.isdir(absin):
            print("The input folder is not correct")
            exit()
        
        docker = os.path.join(absal, 'docker', 'run_docker.py')

        if not os.path.isfile(docker):
            print("There is no run_docker.py file at {}, check the path".format(docker))
            exit()
        else: 
            if debug:
                print("The run_docker.py file is at {}".format(docker))

        # start processing 
        os.chdir(absal)
        if debug:
            print("Current directory is {}".format(os.getcwd()))
        lfiles = os.listdir(absin)
        for file in lfiles:
            if '.fasta' in file:
                if debug:
                    print(os.path.join(absin, file))
                if not os.path.isdir(args['outpath']):
                    print("Creating folder {}".format(args['outpath']))
                    try:
                        os.mkdir(args['outpath'])
                    except OSError as error: 
                        print(error)  
                print("======================Predicting file {}======================".format(os.path.join(absin, file)))
                command = ['python3', 'docker/run_docker.py', ' --output_dir=' + absout, 
                    '--max_template_date=' + args['templatedate'], 
                    '--data_dir=' + args['datadir'], '--run_relax=' + args['relax'],
                    '--model_preset=' + args['modelpreset'],
                    '--fasta_paths=' + os.path.join(absin, file)]
                if debug:
                    print(' '.join(command))
                    print('')
                os.system(' '.join(command))
                
if __name__ == "__main__":
    main()
    