import argparse as ap
from ast import parse
from datetime import date
from logging import exception
import numpy as np 
import os
from datetime import date
import subprocess


def parseAgruments():
    parser = ap.ArgumentParser(description='Superimpose alphafold predicted \
        structure with pdb structure to have coordinate to use with ligands') 
    parser.add_argument('-f', '--folder', type=str, nargs='?', \
        help='Path of the folder contains subfolders, '/ + \
            'each subfolder is a prediction of alphafold for a mutation, ' + \
            'the ranked_0.pdb file is inside the subfolders', \
            required=True)
    parser.add_argument('-o', '--outpath', type=str, nargs='?', \
        help='Path of output directory', required=True)
    parser.add_argument('-c', '--chimerax', type=str, nargs='?', \
         help='Path of ChimeraX running file', required=True) 
    parser.add_argument('-p', '--model', type=str, nargs='?', \
         help='Path of the model pdb file to align', required=True) 
    return parser
            

def main():
    parser = parseAgruments()
    args= vars(parser.parse_args())
    if len(args) < 4:
        parser.print_help()
    else:
        debug = args['debug']
        if debug:
            print('Parsed enough arguments, processing')
        
        # get absolute path of all the input directory 
        absin = os.path.abspath(args['folder'])
        abschi = os.path.abspath(args['chimerax'])
        absout = os.path.abspath(args['outpath']) 
        abspdb = os.path.abspath(args['model'])
        if debug:
            print("The absolute paths:")
            print(absin)
            print(abschi)
            print(absout)
            print(abspdb)
        
        # check if the paths are correct 
        if not os.path.isdir(absin):
            print("The input folder is not correct")
            exit()

        if not os.path.isfile(abschi):
            print("There is no ChimeraX file at {}, check the path".format(abschi))
            exit()
        else: 
            if debug:
                print("The ChimeraX file is at {}".format(abschi))

        # start processing 
        os.chdir('/') # go to root directory of the machine
        if debug:
            print("Current directory is {}".format(os.getcwd()))

        contents = os.scandir(absin)
        for content in contents:
            if os.is_dir (content):
                if debug:
                    print(os.path.join(absin, content))
                
if __name__ == "__main__":
    main()
    