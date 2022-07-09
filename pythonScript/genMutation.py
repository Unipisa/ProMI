import argparse as ap
import numpy as np 
import os


def parseAgruments():
    parser = ap.ArgumentParser(description='Generate mutated sequences') 
    parser.add_argument('-f', '--fasta', type=str, nargs='?', help='Path to fasta file of wild type sequence', required=True)
    parser.add_argument('-m', '--mutations', type=str, nargs='?', help='Path to list of mutations file', required=True)
    parser.add_argument('-o', '--outpath', type=str, nargs='?', help='Path of output directory', required=True)
    parser.add_argument('-d', '--debug', help="Print log to debug or not", action='store_true')
    return parser

def readFasta(path, debug=False):
    if debug:
        print('Process fasta file')

    f=open(path, 'r')
    lines = f.readlines()
    print(lines)
    seqs = []
    des = ''
    seq = ''
    for line in lines:
        if line[0] == '>':
            if debug:
                print('Processing {}'.format(line))
            if seq != '':
                if seq.isalpha():
                    seqs.append({'sequence': seq.upper(), 'description': des})
                    des = ''
                else:
                    print("Sequence {} contains invalid symbol(s)", seq)
                    return []
            des = line.replace('\n', '')
            seq = ''
        else: 
            seq = seq + line.replace('\n', '') 
    # append the last sequence to the list
    if seq != '':
        if seq.isalpha():
            seqs.append({'sequence': seq, 'description': des})
            des = ''
        else:
            print("Sequence {} contains invalid symbol(s)", seq)
            return []
    if debug: 
        print(seqs)
    return seqs 
            

def parseMutationLine(line, debug=False): # format new_aa start_pos stop_pos seq_id 
    words = line.split()
    newaa = words[0]
    start = int(words[1])
    stop = int(words[2])
    if len(words) >= 4:
        seq_id = int (words[3]) - 1
    else:
        seq_id = 0
    if debug:
        print("Mutation {} {} {} {}".format(newaa, start, stop, seq_id + 1))
    return [newaa, start, stop, seq_id]


def parseMutationFile(path, debug=False):
    f = open(path, 'r')
    lines = f.readlines()
    muts = []
    for line in lines:
        mut = parseMutationLine(line, debug)
        muts.append(mut)
    return muts 

def writeNewFile(fastas, fullname, debug=False):
    outfile = open(fullname, 'w')
    for fasta in fastas:
        outfile.write(fasta['description']) 
        outfile.write('\n')
        outfile.write(fasta['sequence'])
        outfile.write('\n')
    outfile.close()


def genMutationFiles(fastas, muts, outpath, debug=False): # start counts from 1
    for mut in muts:
        seq_id = mut[3]
        fasta = fastas[seq_id].copy()
        new_aa = mut[0]
        start = mut[1]
        stop = mut[2]
        if start <= 0:
            print("Start cannot be smaller or equal zero")
            return
        if stop > len(fasta['sequence']):
            print('Stop cannot be larger than the length of the sequence which is {}'.format(len(fasta['sequence'])))
            return
        if stop < start:
            print("Stop should be equal or large than start")
            return 
        for i in range(start-1, stop):
            seq = list(fasta['sequence'])
            des = fasta['description']
            old_aa = seq[i]
            seq[i] = new_aa
            new_fasta = {'sequence': "".join(seq), 'description': des}
            seq = []
            new_fastas = fastas.copy()
            new_fastas[seq_id] = new_fasta 
            name = old_aa + str(i + 1) + new_aa + '.fasta'
            fullname = os.path.join(outpath, name)
            if debug:
                print('Save file to {}'.format(fullname))
            # write new fasta file 
            writeNewFile(new_fastas, fullname, debug)
            

def main():
    parser = parseAgruments()
    args= vars(parser.parse_args())
    if len(args) < 4:
        parser.print_help()
    else:
        debug = args['debug']
        if debug:
            print('Parsed enough arguments, processing')
        seqs = readFasta(args['fasta'], debug)
        if len(seqs) <= 0:
            return
        muts = parseMutationFile(args['mutations'], debug)
        genMutationFiles(seqs, muts, args['outpath'], debug)

    

    



if __name__ == "__main__":
    main()
    