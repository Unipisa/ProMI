import argparse as ap
import os
import calRMSD 
from genShareFolder import parseLine 




def parseAgruments():
    parser = ap.ArgumentParser(description='Superimpose alphafold predicted \
        structure with pdb structure to have coordinate to use with ligands') 
    parser.add_argument('-f', '--folder', type=str, nargs='?', \
        help='Path of the folder contains subfolders, ' + \
            'each subfolder is a prediction of alphafold for a mutation, ' + \
            'the ranked_0.pdb file is inside the subfolders', \
            required=True)
    parser.add_argument('-o', '--outpath', type=str, nargs='?', \
        help='Path of output directory', required=True)
    parser.add_argument('-c', '--chimera', type=str, nargs='?', \
         help='Path of Chimera running file', required=True) 
    parser.add_argument('-p', '--model', type=str, nargs='?', \
         help='Path of the model pdb file to align', required=True) 
    parser.add_argument('-l', '--ligand', type=str, nargs='?', \
                        help='Name of ligand of interest', default=None)
    parser.add_argument('-d', '--debug', help="Print log to debug or not", action='store_true')
    
    return parser
            

def takeOneChain(inpdb, outpdb, debug=False):
    file = open(inpdb, 'r')
    lines = file.readlines() 
    prochains = {} 
    for line in lines:
        if line[0:4] == 'ATOM' or line[0:3] == 'TER' or line[0:6] == "HETATM":
            # parse ATOM line
            linedict = parseLine(line, debug=False)
            chain = linedict['chain']
            if chain not in prochains:
                if debug:
                    print("New protein chain: {}".format(chain))
                prochains[chain] = [linedict['line']]
            else:
                # if debug:
                #     print("Append line to old chain {}".format(chain))
                prochains[chain].append(linedict['line'])
    
    # print(list(prochains.keys())) 
    keylist = list(prochains.keys())
    if len(keylist) >= 1:
        if debug:
            print("Work with only monomer, take only one chain {}".format(keylist[0]))
        
        newfile = open(outpdb,'w')
        newfile.write(''.join(prochains[keylist[0]]))
        newfile.close() 
    else: 
        print("Empt keylist: ")
        print(keylist)

        # for chainline in prochains[keylist[0]]:
        #     # print(chainline['line'])
        #     # print(chainline) 
        #     newfile.writeline

        


def main():
    parser = parseAgruments()
    args= vars(parser.parse_args())
    if len(args) < 4:
        parser.print_help()
    else:
        debug = args['debug']
        if debug:
            print('========================Start processing========================')
            print('Parsed enough arguments, processing')
        
        # get absolute path of all the input directory 
        absin = os.path.abspath(args['folder'])
        abschi = os.path.abspath(args['chimera'])
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
            print("There is no Chimera file at {}, check the path".format(abschi))
            exit()
        else: 
            if debug:
                print("The Chimera file is at {}".format(abschi))

        pdbnameonly = abspdb.split('/')[-1]
        mononame = 'mono_' + pdbnameonly 
        absmononame = abspdb.replace(pdbnameonly, mononame)
        if debug:
            print("Monomer file is saved at {}".format(absmononame))
        takeOneChain(abspdb, absmononame)
        

        # this block is for run chimera from my macbook
        # os.chdir('/') # go to root directory of the machine
        # if debug:
            # print("Current directory is {}".format(os.getcwd()))

        if not os.path.isdir(absout):
            if debug:
                print('Output directory is not exist, create one')
            os.mkdir(absout)
        # os.system('/usr/bin/chimerax -h')
        # start processing 
        contents = os.scandir(absin)
        for content in contents:
            if os.path.isdir(content):
                if debug:
                    print('============================')
                    print(os.path.join(absin, content))
                rank0 = os.path.join(absin, content, 'ranked_0.pdb')
                if os.path.isfile(rank0):
                    if debug:   
                        print('Align the file {} with the file {}'.format(rank0, absmononame))
                    # write command line .cxc file for ChimeraX
                    # store in the subfolder of the mutated structure
                    if not os.path.isdir(os.path.join(str(absout), str(content.name))):
                        if debug:
                            print('Create subfolder {} for output'.format(os.path.join(str(absout), str(content.name))))
                        os.mkdir(os.path.join(str(absout), str(content.name)))
                    
                    # X stuffs are for ChimeraX 
                    filecontentX = 'mmaker #2 to #1\n'
                    filecontentX += 'save ' + str(absout) + '/' + str(content.name)+ '/' + \
                        str(content.name) + '_aligned.pdb ' + \
                        'models #2 relModel #1\n'
                    filecontentX += 'exit'

                    filenameX = os.path.join(absin, content, 'align.cxc')
                    fileX = open(filenameX, 'w')
                    fileX.write(filecontentX)
                    fileX.close()


                    # non X stuffs are for Chimera 
                    filecontent = 'mmaker #1 #0\n'
                    filecontent += 'write format pdb relative #0 #1 ' + \
                        str(absout) + '/' + str(content.name)+ '/' + \
                        str(content.name) + '_aligned.pdb\n'
                    filecontent += 'stop\n'

                    filename = os.path.join(absin, content, 'align.cmd') 
                    file = open(filename, 'w')
                    file.write(filecontent)
                    file.close()

                    if debug:
                        print("Wrote cxc and cmd file")

                    # calling Chimera to align file with the pdb model.
                    # command = ['.' + abschi, '--nogui', '--silent' ,abspdb, str(rank0), filename]

                    command = [abschi, '--nogui', '--silent' ,absmononame, str(rank0), filename]

                    # commandX = ['.' + abschi, '--nogui', '--silent' ,abspdb, str(rank0), filenameX]
                    commandX = [abschi, '--nogui', '--silent' ,absmononame, str(rank0), filenameX]

                    if 'chimerax' in str(abschi).lower():
                        if debug:
                            print(' '.join(commandX))
                            print(filecontentX)
                        os.system(' '.join(commandX))
                    else:
                        if debug:
                            print(' '.join(command))
                            print(filecontent)
                        os.system(' '.join(command))
                    

                    # calculate confidence measurement from here 
                    # between experimental structure (abspdb) 
                    # and aligned alphafold structure (str(content.name) + '_aligned.pdb\n')
                    if args['ligand'] != None:
                        print("Calculate RMSD considering the ligand of interest {}".format(args['ligand']))
                        calRMSD.apply(absmononame, \
                            str(absout) + '/' + str(content.name)+ '/' + str(content.name) + '_aligned.pdb', \
                            str(absout) + '/' + str(content.name)+ '/' + (str(content.name)) + '_confedence.txt', args['ligand'])
                    
                    if debug:
                        print("Done processing, delete {}".format(absmononame))
                    

if __name__ == "__main__":
    main()
    