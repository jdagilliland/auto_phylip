'''
auto_phylip contains functions to handle Phylip programs.
'''
import csv
import os
import subprocess as sub

phy_exec_default = ['phylip','dnapars']
lst_phy_opts_default = ['S','Y','I','4','5','.']

def tab2phy(tabfile, germline=None, outfile=None):
    with open(tabfile,'rb') as f:
        reader = csv.DictReader(f,delimiter='\t')
        lst_dict_entries = [row for row in reader]
    # lst_seqpairs is a list of tuples (sequence id <9-char>, sequence)
    lst_seqpair = [(entry['SEQUENCE_ID'][-9:],entry['SEQUENCE'])
        for entry in lst_dict_entries]
    if germline:
        lst_seqpair = [('Germline', germline)] + lst_seqpair
    n_clone = len(lst_seqpair)
    # WARNING here I assume that all sequences in the tabfile have the same
    # length
    len_clone = len(lst_seqpair[0][1])
    if outfile == None:
        outfile = tabfile.rpartition('.')[0] + '.phy'
    with open(outfile,'wb') as f:
        # write header info
        f.write(_phyrow(n_clone, len_clone))
        # write sequences from tabfile
        for seqpair in lst_seqpair:
            f.write(_phyrow(*seqpair))
    return None

def _phyrow(name, sequence):
    return str(name).ljust(10,' ') + str(sequence) + '\n'

def run_phylip(
    phy_in,
    phy_exec=phy_exec_default,
    lst_phy_opts=lst_phy_opts_default,
    ):
    """
    Run a phylip program with a given phylip executable, and given options.

    The given phylip program is run on the given input file, with the given
    options, and the output files 'outfile' and 'outtree' are renamed to
    `basename`.out, and `basename`.tree respectively.
    `basename` is derived from the input file name.

    Parameters
    ----------
    phy_in : str
        Filename of phylip formatted file to use as input.
    phy_exec : list, optional
        List of strings which are the command line args needed to call
        the desired phylip command, default is ['phylip', 'dnapars']
    lst_phy_opts : list, optional
        List of strings which are the interactive input options for the
        desired phylip program.
        These need to be carefully considered prior to running,
        and must be in the proper order to be meaningful.
        Default is: ['S','Y','I','4','5','.'].
    """
    basename = phy_in.rpartition('.')[0]
    # remove old files
    if os.path.lexists('infile'):
        os.remove('infile')
    if os.path.lexists('outfile'):
        os.remove('outfile')
    if os.path.lexists('outtree'):
        os.remove('outtree')
    # write a command file with the specified options
    cmdfname = '.cmdfile'
    with open(cmdfname,'w') as f:
        f.write(phy_in + '\n')
        for arg in lst_phy_opts:
            f.write(arg + '\n')
        f.write('Y')
    # open phylip process
    p = sub.Popen(phy_exec, stdin=sub.PIPE)
    # send command file contents as input, wait for output
    p.communicate(open(cmdfname,'r').read())
    # rename output files
    try:
        os.rename('outfile', basename + '.out')
        os.rename('outtree', basename + '.tree')
    except:
        print('The expected output was not generated. Phylip may have failed')
        raise
    try:
        os.remove(cmdfname)
    except: raise
    return None
