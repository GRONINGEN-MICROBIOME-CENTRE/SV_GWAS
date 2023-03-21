#!/usr/bin/python
import sys
import argparse
import pandas as pd

"""
Append one or more columns from the second file to the main file matching the files by a column
"""

def parseArguments():
        parser = argparse.ArgumentParser(description="Add column from another file to the current file", add_help=True, epilog = "Finished",conflict_handler="resolve")
        parser.add_argument('-i', type = str, help = "Main file to append columns to", required = True, dest = 'i_fname')
        parser.add_argument('-f', type = str, help = "File from which additional columns are taken", required = True, dest = 'f_fname')
        parser.add_argument('--i_match_col', '-i_m', type = str, help = "Column in the -i file to match files", default = "0", dest = 'i_m')
        parser.add_argument('--f_match_col', '-f_m', type = str, help = "Column in the -f file to match files", default = "0", dest = 'f_m')
        parser.add_argument('--f_add_cols', '-f_cols', type = str, help = "Column in the -f file to match files", required = True, dest = 'f_cols')
        parser.add_argument('--header', action = "store_true", default = False, help = "First lines in both files are headers", dest = 'header')
        #parser.add_argument('-d','-sep', type = str, help = "Field separator", default = "\t", dest = 'sep')
        parser.add_argument('-na','-na', type = str, help = "NA value for lines absent in the f file", default = "NA", dest = 'na')
        parser.add_argument('--col_names','-cn', type = str, help = "new column names separated by commas", default = "", dest = 'cn')

        args = parser.parse_args()
        
        #if args.sep != None:
        #    args.sep = args.sep.decode('string-escape')

        for arg, val in vars(args).items():
            sys.stderr.write(arg + "\t" + str(val) + "\n")
        
        return args

def main(args):
    i_fname = sys.stdin if args['i_fname'] == "stdin" else args['i_fname']
    header_row = 0 if args['header'] else None
    #separator = args['sep']
    separator = "\t"
    
    i_file = pd.read_csv(i_fname, sep = separator, header = header_row)
    f_file = pd.read_csv(args['f_fname'], sep = separator, header = header_row)
    
    if ',' not in args['i_m']:
        i_m = [i_file.columns[int(args['i_m'])]]
        f_m = [f_file.columns[int(args['f_m'])]]
    else:
        i_m = [i_file.columns[int(x)] for x in args['i_m'].split(',')]
        f_m = [f_file.columns[int(x)] for x in args['f_m'].split(',')]
    
    f_cols = [f_file.columns[int(x)] for x in args['f_cols'].split(",")]
    
    if len(args['cn']) > 0:
        colnames = args['cn'].split(',')
    else:    
        colnames = f_cols

    #elif header_row == 0:
    #sys.stderr.write()
    
    sys.stderr.write("i_m = " +  ",".join(map(str, i_m)) + "; f_m = " + ",".join(map(str,f_m)) + "; f_cols = " + ",".join(map(str,f_cols)) + "; colnames = " + ",".join(map(str,colnames)) + "\n")
    
    if len(f_cols) == 1 and len(i_m) == 1:
        if str(colnames[0]) not in i_file.columns:
            newcol = str(colnames[0])
        else:
            newcol = str(colnames[0])+ ".y"
        res_table = i_file.assign(**{newcol : i_file[i_m[0]].map(f_file.set_index(f_m[0])[f_cols[0]])})
    else:
        res_table = i_file.merge(f_file[f_m + f_cols], how='left', left_on = i_m, right_on = f_m)

    res_table.to_csv(sys.stdout, sep="\t", na_rep=args['na'], index = False, header = args['header'])

if __name__ == "__main__":
    args = vars(parseArguments())
    main(args)