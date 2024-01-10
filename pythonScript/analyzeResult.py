import argparse as ap
import csv
import math 


def parseAgruments():
    parser = ap.ArgumentParser(description='Pass result to analyze') 
    parser.add_argument('-a', '--alpha', type=str, nargs='?', \
        help='Path to file contains results of Alphafold structure', required=True)
    parser.add_argument('-e', '--exp', type=str, nargs='?', \
        help='Path to file contains results of experimental structure', required=True)
    parser.add_argument('-d', '--debug', \
        help="Print log to debug or not", action='store_true')
    return parser

def readResFile(path, debug=False):
    alres = {}
    with open(path) as falpha:
        alpha_cont = csv.reader(falpha, delimiter=',')
        for row in alpha_cont:
            if debug:
                print(row)
            name = row[0].strip().lower()
            status = row[1].strip().lower()
            try:
                aff = float(row[2])
            except:
                aff = 0.0
            try:
                var = float(row[3])
            except:
                var = 0.0
            ref = row[4].strip().lower()    

            if name not in alres:
                alres[name] = {} 
                alres[name]['aff'] = []
                alres[name]['var'] = []

            alres[name]['status'] = status
            alres[name]['aff'].append(aff)
            alres[name]['var'].append(var)
            alres[name]['ref'] = ref

    return alres

def calRelVar(res, thres, debug=False): # calculate relative variance and compare with threshold 
    for key, value in res.items():
        # print(key)
        affs = res[key]['aff'].copy()
        mean = sum(affs)/float(len(affs))
        var = 0.0 
        for val in res[key]['aff']:
            var += (val - mean)*(val-mean)
        var = math.sqrt(var)
        var = var/abs(mean) 

        var = var/(len(affs))
        while var > thres: # if relative var greater than threshold, exclude the largest value 
            # if debug:
                # print('Old var is {} and old mean is {}'.format(var, mean))
                # print(key)
                # print("Var {}, threshold {}".format(var, thres))
            maxval = max(affs)
            # print(maxval)
            affs.remove(maxval)
            # print(affs)
            mean = sum(affs)/float(len(affs))
            var = 0.0
            for val in affs:
                var += (val - mean)*(val-mean)
            var = math.sqrt(var)
            var = var/abs(mean) 
            var = var/(len(affs))
            
        res[key]['mean'] = mean
        res[key]['relVar'] = var
        print(key, res[key]['mean'])
    return res 

def getWtMut(res, debug=False):
    listwt = []
    listmut = []
    for key, value in res.items():
        if value['status'] == 'w':
            listwt.append(key)
        else:
            listmut.append(key)
    return listwt, listmut

def checkRes(res, debug=False):
    listcorr = []
    for key, value in res.items():
        # print(key)
        ref = value['ref']
        if ref == 'wt':
            continue
        else:
            wtmean = res[ref]['mean']
            gt = value['status']
            if gt == 'i' and wtmean > value['mean']:
                # print("Mutant mean -  wt mean = {}".format((value['mean'] - wtmean)/wtmean))
                listcorr.append(key)
            elif gt == 'd' and wtmean < value['mean']:
                # print("Mutant mean -  wt mean = {}".format((value['mean'] - wtmean)/wtmean))
                listcorr.append(key)
            # else:
                # print("Wrong: {}".format((value['mean'] - wtmean)/wtmean))
            # print (- wtmean + value['mean'])
    return listcorr 
            

def main():
    parser = parseAgruments()
    args= vars(parser.parse_args())
    if len(args) < 3:
        parser.print_help()
    else:
        debug = args['debug']
        
    alres = readResFile(args['alpha'])
    exres = readResFile(args['exp'])

    listwt, listmut = getWtMut(alres)

    if debug:
        print("List wild type proteins")
        print(listwt)
        print("List mutant proteins")
        print(listmut)

    # for i in (0.05, 0.1, 0.2, 0.3, 0.4, 0.45, 0.5, 0.6, 1.0):
    # for i in [0.05]:
    for i in (1000.0, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.09, 0.08, 0.07, 0.06, 0.05, 0.04, 0.03, 0.02, 0.01, 0):
        print("With threshold {}:".format(i))
        # work with alphafold result 
        print("Alphafold result")
        alres = calRelVar(alres, i, debug)
        alcorr = checkRes(alres, debug)
        print("Number of correct prediction with alphafold {}".format(len(alcorr)))
        print("Accounts for {} percent".format(float(len(alcorr))/float(len(listmut))))

        print("Correct proteins: ")
        print(alcorr)
        alincorr = list(set(listmut) - set(alcorr))
        print("Incorrect proteins: ")
        print(alincorr)


        print('\n')

        print("Experimental result")
        # work with experimental result 
        exres = calRelVar(exres, i, debug)
        excorr = checkRes(exres, debug)
        print("Number of correct prediction with experimental structure {}".format(len(excorr)))
        print("Accounts for {} percent".format(float(len(excorr))/float(len(listmut))))
        
        print("Correct proteins: ")
        print(excorr)
        exincorr = list(set(listmut) - set(excorr))
        print("Incorrect proteins: ")
        print(exincorr)

        print("Number of agreement for correct: {}".format(len(set(alcorr)&set(excorr))))
        print("Number of agreement for incorrect: {}".format(len(set(alincorr)&set(exincorr))))
        print('--------------------')


        print('\n\n')

'''
python3 pythonScript/workflowExpStr.py -w /home/t.pham1/workspace -o scripts/structures_res -i scripts/structures/1w6y_ex/input_1w6y_ff99_15.txt -p scripts/structures/1w6y_ex -f pythonScript/data -c workflow -l acpype/acpype -d
'''   

if __name__ == "__main__":
    main()