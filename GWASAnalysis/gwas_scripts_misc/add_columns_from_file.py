#!/usr/bin/python
"""
Append one or more columns from the second file to the main file matching the files by a column
"""


import sys
import argparse
import gzip

def parseArguments():
        parser = argparse.ArgumentParser(description="Add column from another file to the current file", add_help=True, epilog = "Finished",conflict_handler="resolve")
        parser.add_argument('-i', type = str, help = "Main file to append columns to", required = True, dest = 'i_fname')
        parser.add_argument('-f', type = str, help = "File from which additional columns are taken", required = True, dest = 'f_fname')
        parser.add_argument('--i_match_col', '-i_m', type = int, help = "Column in the -i file to match files", default = 0, dest = 'i_m')
        parser.add_argument('--f_match_col', '-f_m', type = int, help = "Column in the -f file to match files", default = 0, dest = 'f_m')
        parser.add_argument('--f_add_cols', '-f_cols', type = str, help = "Column in the -f file to match files", required = True, dest = 'f_cols')
        parser.add_argument('--header', action = "store_true", default = False, help = "First lines in both files are headers", dest = 'header')
        parser.add_argument('-d','-sep', type = str, help = "Field separator", default = "\t", dest = 'sep')
        parser.add_argument('-na','-na', type = str, help = "NA value for lines absent in the f file", default = "NA", dest = 'na')

        args = parser.parse_args()
        
        if args.sep != None:
            args.sep = args.sep.decode('string-escape')

        for arg, val in vars(args).items():
            sys.stderr.write(arg + "\t" + str(val) + "\n")
        
        return args

def main(args):
    sep = args['sep']
    f_cols = args['f_cols'].split(",")
    i_m = args['i_m']
    f_m = args['f_m']

    add_dict = {}
    add_col_names = []
    
    if f_cols == ["all"]:
        f_line = open(args['f_fname']).readline().rstrip("\r\n").split(sep)
        ncols = len(f_line)
        f_cols = range(1,ncols)
    with open(args['f_fname']) as f_file:
        if args['header']:
            header = f_file.readline().rstrip("\r\n").split(sep)
            add_col_names = [header[int(i)] for i in f_cols]
        for line in f_file:
            spl = line.rstrip("\r\n").split(sep)
            if spl[f_m] in add_dict:
                sys.stderr.write ("Duplicate row names in the second file are not allowed: " + spl[f_m] + "\nExiting!")
                sys.exit(-1)
            add_dict[spl[f_m]] = [spl[int(i)] for i in f_cols]


    if args['i_fname'] != "stdin":
        if args['i_fname'].endswith(".gz"):
            i_file = gzip.open(args['i_fname'])
        else:
            i_file = open(args['i_fname'])
    else:
        i_file = sys.stdin
    
    if args['header']:
            print (i_file.readline().rstrip("\r\n") + sep + sep.join(add_col_names))
    
    for line in i_file:
        spl = line.rstrip("\r\n").split(sep)
        add_col_list = add_dict.get(spl[i_m], [args['na']] * len(f_cols))
        print (line.rstrip("\r\n") + sep + sep.join(add_col_list))



if __name__ == "__main__":
    args = vars(parseArguments())
    main(args)