'''
auto_phylip contains functions to handle Phylip programs.
'''
import csv
import os
import subprocess as sub
import re

phy_exec_default = ['phylip','dnapars']
lst_phy_opts_default = ['S','Y','I','4','5','.']

def tab2phy(lst_tabfile, germline=None, outfile=None, **kwarg):
    match = kwarg.pop('match', None)
    column = kwarg.pop('column', 'CLONE')
    flags = kwarg.pop('flags', 0)
    lst_dict_entries = _gather_entries(lst_tabfile)
    lst_dict_entries = _filter_entries(lst_dict_entries,
        match=match, column=column, flags=flags)
    if outfile == None and len(lst_tabfile) == 1:
        outfile = lst_tabfile[0].rpartition('.')[0] + '.phy'
    elif outfile == None:
        outfile = 'file.phy'
    lst_entries2phy(lst_dict_entries, outfile)
    return None

def _entry2seqpair(entry):
    return (entry['SEQUENCE_ID'][-9:],entry['SEQUENCE'])

def lst_entries2phy(lst_dict_entries, outfile, **kwarg):
    """
    Write a list of tabfile entries `lst_dict_entries` to a PHYLIP formatted
    `outfile`.

    Parameters
    ----------
    lst_dict_entries : list
        A list of dictionaries which each represent a row of data from a
        tabfile.
    outfile : str
        A filename to which to output the PHYLIP formatted data.
        This is mandatory at this level, since there is no other way to
        handle coming up with an output filename given the args.
    germline : str, optional, no-op

    Returns
    -------
    None
    """
    germline = kwarg.pop('germline', None)
    # lst_seqpairs is a list of tuples (sequence id <9-char>, sequence)
    lst_seqpair = [_entry2seqpair(entry) for entry in lst_dict_entries]
    if germline:
        lst_seqpair = [('Germline', germline)] + lst_seqpair
    else:
        dict_germline_seq = dict(
            [(entry['CLONE'], entry.get('GERMLINE_GAP_DMASK'))
            for entry in lst_dict_entries]
            )
    n_seq = len(lst_seqpair)
    if n_seq == 0:
        raise ValueError('''No matches found.''')
    len_seq = len(lst_seqpair[0][1])
    for seqpair in lst_seqpair:
        if len(seqpair[1]) != len_seq:
            raise ValueError('''Not all sequences are of the same length.''')
    with open(outfile,'wb') as f:
        # write header info
        f.write(_phyrow(n_seq, len_seq))
        # write sequences from tabfile
        for seqpair in lst_seqpair:
            f.write(_phyrow(*seqpair))
    return None
def _gather_entries(lst_file):
    """
    Gather tabfile entries from a list of tabfiles.
    """
    lst_entries = list()
    for tabfile in lst_file:
        with open(tabfile,'rb') as f:
            reader = csv.DictReader(f,delimiter='\t')
            lst_entries.extend([row for row in reader])
    return lst_entries

def _filter_entries(lst_dict_entries, match=None, column='CLONE', flags=0):
    """
    Filter a list of tab file entries `lst_dict_entries` by matching
    the text in a give column `column` to the regex string provided in
    `match`.

    Parameters
    ----------
    lst_dict_entries : list
        A list of dictionaries which each represent a row of data from a
        tabfile.
    match : str
        A regex string to match entries.
    column : str, optional
        The column whose data to match (default: 'CLONE').
    flags : int, optional
        An int representing the flags to apply for the regex matching
        (default: 0, no flags).

    Returns
    -------
    lst_dict_entries_match : list
        A list of dict entries filtered by the selected column matching the
        regex.
    """
    if type(match) == str:
        reg_match = re.compile(match, flags)
#    elif type(match) == _sre.SRE_Pattern:
#        # check if the provided match is a compiled regex, if so, use as is
#        reg_match = match
    elif match == None:
        # if match not provided, return the supplied lst_dict_entries
        # unchanged
        return lst_dict_entries
    else:
        print(match)
        print(type(match))
        raise TypeError('''`match` is of inappropriate type, must be `str`
            or compiled regex.''')
    # I use the match function, because I want users to be specific about
    # whether they want to match the beginning of the string or w/e.
    lst_dict_entries_match = [entry for entry in lst_dict_entries if
        reg_match.match(entry.get(column))]
    return lst_dict_entries_match

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
