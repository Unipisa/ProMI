import argparse as ap
import csv
import math 
import os
from genShareFolder import parseLine 


def parseAgruments():
    parser = ap.ArgumentParser(description='Pass result to analyze') 
    parser.add_argument('-p', '--pdb', type=str, nargs='?', \
        help='Path to pdb file', required=True)
    parser.add_argument('-o', '--output', type=str, nargs='?', \
        help='Path to folder for the output file', required=True)
    parser.add_argument('-d', '--debug', \
        help="Print log to debug or not", action='store_true')
    return parser


def callPropka(pdb, debug=False):
    if not os.path.isfile(pdb):
        print("There is no pdb file to run at {}".format(pdb))
        exit()
    call = ['propka3', pdb, '-o', '7.4']
    os.system(' '.join(call)) 
    pwd = os.path.abspath(os.getcwd())
    outname = os.path.join(pwd, pdb.split('/')[-1].replace('.pdb', '.pka'))
    print("Generate a .pka file at {}".format(outname))
    return outname 
    
    # output will be generated at the current folder 
    # name of output will be name of input file with suffix '.pka'



def parsePkafile(pka, debug=False):
    aas =["HIS", "ASP", "GLU", "CYS", "LYS", "ARG", "TYR"]
    if not os.path.isfile(pka):
        print("Cannot find pka file")
        exit()
    fpka = open(pka, 'r')
    lines = fpka.readlines()
    start = False
    stop = False 
    reslines = []
    for line in lines:
        if "SUMMARY OF THIS PREDICTION" in line:
            start = True 
        if start and not stop:
            if '------------------------------------------------' in line:
                stop = True
            else:
                reslines.append(line)
        if stop:
            break 
    
    aminos = []
    for line in reslines:
        # if debug:
            # print(line)
        words = line.split()
        for aa in aas:
            if aa in words and len(words) == 5:
                # print(aa, words)
                amino = dict()
                amino['name'] = words[0]
                amino['id'] = words[1]
                amino['chain'] = words[2]
                amino['pka'] = words[3]
                amino['mpka'] = words[4]
                aminos.append(amino)
                # print(amino)
    # print(aminos)
    # for amino in aminos:
    #     print(amino)
        #delete propka file 
    if debug:
        print("Done processing, remove file {}".format(pka))
    os.remove(pka)

    return aminos
        
def modifyPdb(pdb, aminos, ph=7.4, debug=False):
    paa = { "HIS":{'deprot':"HIE", "prot": "HIS", 'dedeprot':"HID"}, \
            "ASP":{'deprot':"ASP", 'prot': "ASH"}, \
            "GLU":{'deprot':"GLU", 'prot': "GLH"}, \
            "CYS":{'deprot':"CYM", 'prot': "CYS"}, \
            "LYS":{'deprot':"LYN", 'prot': "LYS"}, \
            "TYR":{'deprot':"TYR", 'prot': "TYR"}, \
            "ARG":{'deprot':"ARG", 'prot': "ARG"}}
    for amino in aminos:
        # case of HIS 
        if amino['name'] == 'HIS':
            if float(amino['pka']) < 4.0:
                amino['name'] = paa['HIS']['dedeprot']
            elif float(amino['pka']) >= 4.0 and float(amino['pka']) < ph:
                amino['name'] = paa['HIS']['deprot']
            else:
                amino['name'] = paa['HIS']['prot']
        else:
            if float(amino['pka']) > ph: # use protonated 
                amino['name'] = paa[amino['name']]['prot']
            else:
                amino['name'] = paa[amino['name']]['deprot']
            # print(paa[amino['name']])

    fpdb = open(pdb, 'r')
    lines = fpdb.readlines()
    linedicts = []

    protfname = 'prot_' + pdb.split('/')[-1]
    
    fullprotname = pdb.replace(pdb.split('/')[-1], protfname)
    fout = open(fullprotname, 'w')
                          
    for line in lines:
        # print(line)
        if line[0:4] == 'ATOM' or line[0:3] == 'TER':
            linedict = parseLine(line)
            for amino in aminos:
                if int(linedict['rid']) == int(amino['id']):
                    linedict['line'] = linedict['line'].replace(str(linedict['name']), str(amino['name']))
                    linedict['name'] = amino['name']
                    # print(linedict['line'])
            linedicts.append(linedict)
            fout.write(linedict['line'])
        else:   
            fout.write(line)

    fout.close()
    if debug:
        print("Save protonated file at {}".format(fullprotname))
    

def main():
    parser = parseAgruments()
    args= vars(parser.parse_args())
    if len(args) < 3:
        parser.print_help()
    else:
        debug = args['debug']
    pka = callPropka(args['pdb'], debug)
    aminos = parsePkafile(pka, debug)
    modifyPdb(args['pdb'], aminos, ph=7.4,debug=debug)

if __name__ == "__main__":
    main()