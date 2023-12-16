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
    parser.add_argument('-d', '--debug', help="Print log to debug or not", action='store_true')
    return parser



def main():
    parser = parseAgruments()
    args= vars(parser.parse_args())
    if len(args) < 2:
        parser.print_help()
    else:
        debug = args['debug']
    
    rootpath = args['folder']
    # subfolders = os.scandir(rootpath)
    start = [4900]#, 
    # start = [200, 300, 400, 500, 600, 700, 800, 2000, \
    #     2300, 2600, 2900, 3200, 3500, 4100, 4500, 4700, 4900]#, 5000, \
        # 5400, 5800, 6200, 6600, 6900, 7300, 7600, 7700, 7800, \
        # 8300, 8700, 9100, 9400, 9900, 10300, 10800, 11100, 11500, \
        # 12100, 12500, 12900, 13300, 13800, 14200, 14600, 15000, \
        # 15700, 16100, 16700, 17100, 17700, 18200, 18800, 19000, \
        # 19100, 19200, 19300, 19400, 19500, 19600, 19700, 19800, \
        # 19900]
    stop  = [5000]#, 
    # stop = [300, 400, 500, 600, 700, 800, 900, 2100, \
    #     2400, 2700, 3000, 3300, 3600, 4200, 4600, 4800, 5000]#, 5100, \
        # 5500, 5900, 6300, 6700, 7000, 7400, 7700, 7800, 7900, \
        # 8400, 8800, 9200, 9500, 10000, 10400, 10900, 11200, 11600, \
        # 12200, 12600, 13000, 13400, 13900, 14300, 14700, 15100, \
        # 15800, 16200, 16800, 17200, 17800, 18300, 18900, 19100, \
        # 19200, 19300, 19400, 19500, 19600, 19700, 19800, 19900, \
        # 20000]
    if len(start) != len(stop):
        print("Start and stop are not the same size")
        exit()
    for i in range(len(start)):
        # if i == 0:
        #     continue
        # subfolders = os.scandir(rootpath)
        if debug:
            print("Range {} - {}".format(start[i], stop[i]))
        # modify mmpbsa file inside Share folder 
        if not os.path.isdir(os.path.join(rootpath, 'Share')):
            print("There is no Share folder at {}".format(os.path.join(rootpath)))
        mmname = os.path.join(rootpath, 'Share', 'mmpbsa.in')
        # if not os.path.isfile(mmname):
        #     print("There is not mmpbsa.in file at {}".format(mmname))
        mmfile = open(mmname, 'w')
        mmfile.write('# General namelist variables\n')
        mmfile.write('&general\n')
        mmfile.write('sys_name             = "Prot-Lig-ST" \n')
        mmfile.write('startframe           = ' + str(start[i]) + ' \n')
        mmfile.write('endframe             = ' + str(stop[i]) + ' \n')
        mmfile.write('interval             = 5 \n')
        mmfile.write('verbose              = 1 \n')
        mmfile.write('/\n')

        mmfile.write('&gb\n')
        mmfile.write('igb                  = 2 \n')
        mmfile.write('saltcon              = 0.100\n')
        # mmfile.write('ifqnt                 = 1\n')
        # mmfile.write('qm_residues           = "within 5" \n')
        # mmfile.write('qm_theory             = "PM3"\n')
        mmfile.write('/\n')
        
        mmfile.write('&decomp\n')
        mmfile.write('idecomp               = 2\n')
        mmfile.write('dec_verbose           = 3\n')
        mmfile.write('print_res             = "within 4"\n')
        mmfile.write('/\n')


        mmfile.close()
        com = ['docker', 'run', '--gpus', 'all', \
            '--rm', '--user', '$(id -u):$(id -g)', '--mount',
            'src=/home/t.pham1/workspace,target=/tmp/workspace,type=bind', \
            'workflow', '/tmp/workspace/scripts/bashScript/run_MDBatch.sh']
        os.system(' '.join(com))

        curdir = os.getcwd()
        print ("Current dir {}".format(curdir))

        for subfolder in ['Y190F', 'Y190Y']:
            print (os.path.join(rootpath, subfolder, 'FINAL_RESULTS_MMPBSA.dat'))
            newname = 'FINAL_RESULTS_MMPBSA' + '_' + str(start[i]) + '-' + str(stop[i]) + '.dat'
            com2 = ['mv', os.path.join(rootpath, subfolder, 'FINAL_RESULTS_MMPBSA.dat'),\
                    os.path.join(rootpath, subfolder, newname)]
            os.system(' '.join(com2))


if __name__ == "__main__":
    main()