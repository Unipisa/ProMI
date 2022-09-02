import argparse as ap
import numpy as np 
import os
import re


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
        

def parseRangeMutation(line, debug=False): # format new_aa start_pos stop_pos seq_id 
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

def writeNewFile(fastas, fullname, debug=False):
    outfile = open(fullname, 'w')
    for fasta in fastas:
        outfile.write(fasta['description']) 
        outfile.write('\n')
        outfile.write(fasta['sequence'])
        outfile.write('\n')
    outfile.close()

def genPointMutation(fastas, outpath, line, debug=False): # each line is one mutated fasta file. 
    changes = line.split(';')
    name = ''
    new_fastas = fastas.copy()

    for change in changes: # each mutation is in the form like 1_A45B 
        change = change.strip()
        parts = change.split('_')
        seqid = int(parts[0]) - 1 # the id in the input file starts at 1 while the index in code starts at 0 
        
        if seqid >= len(new_fastas):
            print("There are only {} sequences, \
                meanwhile the sequence index is {}".format(len(new_fastas), seqid))

        fasta = new_fastas[seqid].copy() # fasta is now a dictionary with 2 keys of one sequence 
        seq = list(fasta['sequence']) # seq is now list of characters  

        mut = parts[1]
        pos = re.findall(r'\d+', mut)
        if len(pos) != 1:
            print('Wrong mutation {}'.format(mut))
        pos = int(pos[0])
        aas = mut.split(str(pos))
        pos  -= 1 
        src = aas[0] 
        des = aas[1] 
        if debug:
            print("Consider the mutations: {} {} {}".format(src, pos + 1, des)) 

        if seq[pos] != src:
            print("The original amino acid at position \
                {} is not {} but {}".format(pos + 1, src, seq[pos]))

        seq[pos] = des # change the amino acid here 
        seq = ''.join(seq) # now seq is a string 
        new_des = fasta['description']
        new_des += str(parts[1])
        name += str(parts[1])
        new_fastas[seqid] = {'sequence': seq, 'description':  new_des}

    fullname = os.path.join(outpath, name + '.fasta')
    writeNewFile(new_fastas, fullname, debug)

def parseMutationFile(path, fastas, outpath, debug=False):
    f = open(path, 'r')
    lines = f.readlines()
    muts = []
    for line in lines:
        if 'range' in line:
            line = line.split(':')[1].strip()
            mut = parseRangeMutation(line, debug)
            muts.append(mut)
        if 'point' in line:
            line = line.split(':')[1].strip()
            genPointMutation(fastas, outpath, line, debug)
    return muts 

def writeNewFile(fastas, fullname, debug=False):
    outfile = open(fullname, 'w')
    for fasta in fastas:
        outfile.write(fasta['description']) 
        outfile.write('\n')
        outfile.write(fasta['sequence'])
        outfile.write('\n')
    outfile.close()

def genRangeMutations(fastas, muts, outpath, debug=False): # start counts from 1
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
        muts = parseMutationFile(args['mutations'], seqs, args['outpath'], debug)
        genRangeMutations(seqs, muts, args['outpath'], debug)

if __name__ == "__main__":
    main()