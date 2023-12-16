import os 

infolder = '/Users/phamgiang/Documents/Study/Master.nosync/Thesis/Workspace.nosync/scripts/exp_data/nofix2/'
outfolder = '/Users/phamgiang/Documents/Study/Master.nosync/Thesis/Workspace.nosync/scripts/exp_data/fix/'

files = os.scandir(infolder)
for file in files:
    if '.pdb' in file.name:
        filename = os.path.join(infolder, file.name)
        infile = open(filename, 'r')
        inlines = infile.readlines()

        temname = os.path.join(outfolder, 'fix_' + file.name)
        temfile = open(temname, 'r')
        temlines = temfile.readlines()

        noend = temlines[0:len(temlines) - 1] # omit the last line "END"

        for line in inlines:
            if line[0:6] == 'HETATM':
                noend.append(line)
        
        print(len(temlines))
        noend.append(temlines[len(temlines) - 1])

        outfile = open(os.path.join(outfolder, "fix_" + file.name), 'w')
        outfile.write(''.join(noend))
        outfile.close()
