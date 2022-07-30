import pickle
import argparse as ap
import os



def parseAgruments():
    parser = ap.ArgumentParser(description='Query PAE score') 
    parser.add_argument('-f', '--pkl', type=str, nargs='?', \
        help='Path to .pkl file', required=True)
    parser.add_argument('-x', '--x', type=int, nargs='?', \
        help='X position')
    parser.add_argument('-y', '--y', type=int, nargs='?', \
        help='Y position')
    parser.add_argument('-d', '--debug', \
        help="Print log to debug or not", action='store_true')
    return parser

def main():
    parser = parseAgruments()
    args= vars(parser.parse_args())
    if len(args) < 4:
        parser.print_help()
    else:
        debug = args['debug']
    if not os.path.isfile(args['pkl']) and '.pkl' not in args['pkl']:
        print("File path is not correct")
        exit()

    with open(args['pkl'], 'rb') as f:
        data = pickle.load(f)
        paes = data['predicted_aligned_error']
        # print(paes)
        x = args['x']
        y = args['y']
        try:
            print(paes[x-1][y-1])
        except:
            print("The queried position is out of range") 


if __name__ == "__main__":
    main()


