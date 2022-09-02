import argparse as ap
from ast import parse
from datetime import date
from hashlib import new
from logging import exception
import numpy as np 
import os
from datetime import date
import subprocess


def parseAgruments():
    parser = ap.ArgumentParser(description='Preprocess ligands for MD simulation') 
    parser.add_argument('-f', '--folder', type=str, nargs='?', \
        help='Path to .acpype folder, which is the result of acpype', \
            required=True)
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
            print('========================Start processing========================')
            print('Parsed enough arguments, processing')
        
        # get absolute path of all the input directory 
        absin = os.path.abspath(args['folder'])
        if debug:
            print("The working directory is: {}".format(absin))
        # check if the paths are correct 
        if os.path.isdir(absin):
            if '.acpype' not in str(absin):
                print("The input folder is not correct")
                exit()

        # start processing 
        contents = os.listdir(absin)
        for content in contents:
            # fix topology and parameter files
            if os.path.isfile(os.path.join(absin,content)) and '_GMX.itp' in content:
                if debug:
                    print("Fixing topology and parameter files")
                ligname = content.replace('_GMX.itp','')
                file = open(os.path.join(absin, content))
                lines = file.readlines()
                start, stop = 0, 0
                counter = 0
                for line in lines:
                    if 'atomtypes' in line:
                        start = counter 
                    if 'moleculetype' in line: 
                        stop = counter 
                    counter += 1
                if debug:
                    print ("Ligand's name is: {}".format(ligname))
                prmlines = lines[start:stop]
                iptlines = lines[stop:]
                prmname = os.path.join(absin, ligname + '_fix.prm')
                if debug:
                    print("Save .prm file at: {}".format(prmname))
                prmfile = open(prmname, 'w')
                prmfile.write(''.join(prmlines))

                itpname = os.path.join(absin, ligname + '_fix.itp')
                if debug:
                    print("Save .itp file at: {}".format(itpname))
                itpfile = open(itpname, 'w')
                itpfile.write(''.join(iptlines))

            #fix gromacs coordinate file 
            if os.path.isfile(os.path.join(absin, content)) and '_GMX.gro' in content:
                if debug:
                    print("Fixing .gro file")
                ligname = content.replace('_GMX.gro','')
                file = open(os.path.join(absin, content))
                lines = file.readlines()
                lines_to_fix = lines[2:len(lines) - 1]
                # print(lines_to_fix)
                lines_to_save = [] 
                lines_to_save.append(lines[0])
                lines_to_save.append(lines[1])
                for line in lines_to_fix:
                    words = line.split()
                    new_line = ''
                    for i in range(8 - (len(words[0]) + len(words[1]))):
                        new_line += ' '
                    new_line += words[0] + words[1] # the number of blankspace is important
                    for i in range(7-len(words[2])):
                        new_line += ' '
                    new_line += words[2]
                    for i in range(5 - len(words[3])):
                        new_line += ' '
                    new_line += words[3]
                    for i in range(8 - len(words[4])):
                        new_line += ' '
                    new_line += words[4]
                    for i in range(8 - len(words[5])):
                        new_line += ' '
                    new_line += words[5]
                    for i in range(8 - len(words[6])):
                        new_line += ' '
                    new_line += words[6]
                    new_line += '\n'
                    lines_to_save.append(new_line)
                # print(lines_to_save)
                lines_to_save.append(lines[len(lines)-1])
                groname = os.path.join(absin, ligname + '_fix.gro')
                if debug:
                    print("Save .prm file at: {}".format(groname))
                grofile = open(groname, 'w')
                grofile.write(''.join(lines_to_save))

if __name__ == "__main__":
    main()
    