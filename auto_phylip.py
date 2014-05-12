'''
auto_phylip contains functions to handle Phylip programs.
'''
import csv
import os
import subprocess as sub
import re

phy_exec_default = ['phylip','dnapars']
boot_exec_default = ['phylip', 'seqboot']
cons_exec = ['phylip', 'consense']

n_bootstrap_default = 1000

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
    print('Found {n_match} matches for {str_match}'.format(
        n_match=len(lst_dict_entries), str_match=match))
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
        lst_entries.extend(_get_entries(tabfile))
    return lst_entries

def _get_entries(tabfile):
    """
    Get tabfile entries from a tabfile
    """
    with open(tabfile,'rb') as f:
        reader = csv.DictReader(f,delimiter='\t')
        lst_entries = [row for row in reader]
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
    **kwarg):
    """
    Run a phylip program with a given phylogeny generating executable,
    and given options.

    The given phylogeny program is run on the given input file, with the given
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
    bootstrap : int, optional
        The number of bootstrap iterations for seqboot.
        If not specified, seqboot is skipped.
    """
    # if phy_exec is provided, or None, use default
    phy_exec = kwarg.pop('phy_exec', None)
    if phy_exec == None:
        phy_exec = phy_exec_default
    # it bootstrap is provided use, otherwise None
    bootstrap = kwarg.pop('bootstrap', None)
    if bootstrap:
        phy_in = run_seqboot(phy_in, bootstrap)
    basename = phy_in.rpartition('.')[0]
    # remove old files
    if os.path.lexists('infile'):
        os.remove('infile')
    if os.path.lexists('outfile'):
        os.remove('outfile')
    if os.path.lexists('outtree'):
        os.remove('outtree')
    # write a command file with the specified options
    lst_phy_opts = _get_phy_opts(phy_in, bootstrap)
    cmdfname = write_cmdfile(lst_phy_opts, trailing_nl=False)
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
    if bootstrap:
        # run consense only if bootstrapping was performed
        run_consense(basename + '.tree', **kwarg)
    return None

def run_seqboot(fname, n_bootstrap, **kwarg):
    """
    Run seqboot on a given file for a given number of bootstraps.
    """
    boot_exec = kwarg.pop('boot_exec', boot_exec_default)
    seqboot_opts = _get_seqboot_opts(fname, n_bootstrap, **kwarg)
    # remove old files
    if os.path.lexists('outfile'):
        os.remove('outfile')
    # write a command file with the specified options
    cmdfname = write_cmdfile(seqboot_opts)
    # open phylip process
    p = sub.Popen(boot_exec, stdin=sub.PIPE)
    # send command file contents as input, wait for output
    p.communicate(open(cmdfname,'r').read())
    basename = fname.rpartition('.')[0]
    bootname = basename + '.boot.phy'
    try:
        os.rename('outfile', bootname)
    except:
        print('The expected output was not generated. Phylip may have failed')
        raise
    try:
        os.remove(cmdfname)
    except: raise
    return bootname

def run_consense(fname, **kwarg):
    """
    Run consense on a given set of bootstrapped trees to form a consensus tree.
    """
    consense_opts = _get_consense_opts(fname, **kwarg)
    # remove old files
    if os.path.lexists('intree'):
        os.remove('intree')
    if os.path.lexists('outfile'):
        os.remove('outfile')
    if os.path.lexists('outtree'):
        os.remove('outtree')
    # write a command file with the specified options
    cmdfname = write_cmdfile(consense_opts)
    # open phylip process
    p = sub.Popen(cons_exec, stdin=sub.PIPE)
    # send command file contents as input, wait for output
    p.communicate(open(cmdfname,'r').read())
    basename = fname.rpartition('.')[0]
    try:
        os.rename('outfile', basename + '.cons.out')
        os.rename('outtree', basename + '.cons.tree')
    except:
        print('The expected output was not generated. Phylip may have failed')
        raise
    try:
        os.remove(cmdfname)
    except: raise
    return None

def write_cmdfile(opts, trailing_nl=False):
    """
    Writes a command file for use with PHYLIP programs based on the supplied
    list of options.

    Parameters
    ----------
    opts : list
        A list of strings to write to the commandfile separated by newlines.
    trailing_nl : bool
        Whether or not to add a trailing newline character at the end of the
        file (default: False).

    Notes
    -----
    The last option in `opts` will be written without a trailing newline
    unless otherwise specified.
    This is appropriate behavior for _most_ PHYLIP programs.
    """
    cmdfname = '.cmdfile'
    with open(cmdfname,'w') as f:
        for arg in opts[:-1]:
            f.write(arg + '\n')
        f.write(opts[-1])
    return cmdfname

def _get_phy_opts(fname, bootstrap, **kwarg):
    # SLOPPY: hardcoded arbitrary constant
    seed = 7 # must be odd for some reason!?
    # SLOPPY: hardcoded arbitrary constant
    jumble = 10
    opts = list()
    # first the filename
    opts.append(fname)
    # set search option
    #-----------------
    opts.append('S')#|
    opts.append('Y')#|
    #-----------------
    # bootstrapped datasets will always be interleaved (AFAIK)
    if not bootstrap:
        # specify sequential rather than interleaved data
        opts.append('I')
    if bootstrap:
        # specify multiple datasets
        opts.append('M')
        # specify that the data is multiple, not the weights (right now I
        # don't even use weights)
        opts.append('D')
        # specify the number of datasets, the same as the number of
        # bootstraps
        opts.append(str(bootstrap))
        # provide the seed
        opts.append(str(seed))
        opts.append(str(jumble))
    # these options make the output more rich and thorough
    #-----------------
    opts.append('4')#|
    opts.append('5')#|
    opts.append('.')#|
    #-----------------
    # confirm options
    opts.append('Y')
    return opts

def _get_seqboot_opts(fname, n_bootstrap, weights=False, seed=7):
    lst_seqboot_opts = []
    lst_seqboot_opts.append(fname)
    lst_seqboot_opts.append('R')
    lst_seqboot_opts.append(str(n_bootstrap))
    # the following to be properly implemented later
    if weights:
        lst_seqboot_opts.append('W')
        pass
    # finish putting in options, and confirm
    lst_seqboot_opts.append('Y')
    lst_seqboot_opts.append(str(seed))
    return lst_seqboot_opts

def _get_consense_opts(fname, **kwarg):
    cons_type = kwarg.pop('type', 'mre')
    opts = list()
    opts.append(fname)
    # select consensus type
    if cons_type == 'mre':
        pass
    elif cons_type == 'strict':
        # not implemented
        pass
    elif cons_type == 'mr':
        # not implemented
        pass
    elif cons_type == 'ml':
        # not implemented
        pass
    else:
        print(cons_type)
        raise ValueError(
                'Invalid consensus type {cons}'.format(cons=cons_type))
    # confirm options
    opts.append('Y')
    return opts

def _run_phylip_main():
    import argparse
    parser = argparse.ArgumentParser(
        description=("Run PHYLIP program to build trees from sequences in a" +
            " PHYLIP formatted file."),
        )
    parser.add_argument('-c', '--command',
            dest='command',
            default=None,
            help="""
            Input an string (quoted for 2 or more space separated args),
            which will be the base PHYLIP executable for the runner.
            e.g. 'phylip dnapars' to use DNA parsimony (the most tested
            by the author of this runner),
            or 'phylip dnacomp' to use DNA compatibility,
            or 'phylip dnapenny' to use branch and bound search for trees
            (note, dnapenny may not be practical for more than 10 or 11
            sequences),
            'phylip dnaml' to use DNA maximum likelihood.
            Not all PHYLIP programs that one might choose will output trees,
            or exhibit the predictable behavior on which this script depends.
            """
            )
    parser.add_argument('-b', '--bootstrap',
            dest='bootstrap',
            const=n_bootstrap_default, # default number of replicates
            nargs='?',
            default=None,
            type=int,
            help="""
            Provide a number of bootstrap replicates to use for phylip
            seqboot.
            If you do not include this flag, seqboot will not be used (i.e.
            1 bootstrap).
            If you include this flag with an integer argument,
            that number of bootstrap replicates will be used.
            Be careful with large number of bootstraps, since the phylogeny
            inference will have to be run that much more.
            If you include this flag without any argument, a default value
            will be used.
            (default is {default})
            """.format(default=n_bootstrap_default),
            )
    parser.add_argument('files', nargs='+')
    argspace = parser.parse_args()
    if argspace.command == None:
        lst_cmd_arg = None
    else:
        lst_cmd_arg = argspace.command.split(' ')
    print('''Using {args} to run PHYLIP'''.format(args=lst_cmd_arg))
    for fname in argspace.files:
        run_phylip(fname,
                phy_exec=lst_cmd_arg,
                bootstrap=argspace.bootstrap,
                )
def _run_seqboot_main():
    import argparse
    parser = argparse.ArgumentParser(
            description="""Bootstrap a set of sequence data using phylip
            seqboot""",
            )
    parser.add_argument('-b', '--bootstrap',
            dest='bootstrap',
            const=n_bootstrap_default,
            nargs='?',
            default=None,
            type=int,
            help="""
            Provide a number of bootstrap replicates to use for phylip
            seqboot.
            If you do not include this flag, seqboot will not be used (i.e.
            1 bootstrap).
            If you include this flag with an integer argument,
            that number of bootstrap replicates will be used.
            Be careful with large number of bootstraps, since the phylogeny
            inference will have to be run that much more.
            If you include this flag without any argument, a default value
            will be used.
            (default is {default})
            """.format(default=n_bootstrap_default),
            )
    parser.add_argument('files', nargs='+')
    argspace = parser.parse_args()
    for fname in argspace.files:
        run_seqboot(fname, argspace.bootstrap,
                )
    return None

def _run_consense_main():
    import argparse
    parser = argparse.ArgumentParser(
            description='''Build a consensus tree from bootstrapped trees''',
            )
    parser.add_argument('files', nargs='+')
    argspace = parser.parse_args()
    for fname in argspace.files:
        run_consense(
                fname,
                )
    ### resume here
    return None

def _tab2phy_main():
    import argparse
    parser = argparse.ArgumentParser(
        description='Build a phylip formatted file from a tabfile.',
        )
    parser.add_argument('files', nargs='+')
    parser.add_argument('-m','--match',
        dest='match',
        default=None,
        )
    parser.add_argument('-f','--field',
        dest='column',
        default='CLONE',
        )
    parser.add_argument('-c','--combine',
        dest='phyfname',
        nargs='?',
        const='file.phy',
        default=None,
#        default='file.phy',
        )
    argspace = parser.parse_args()
    tab2phy(
        argspace.files,
        outfile=argspace.phyfname,
        match=argspace.match,
        column=argspace.column,
        )
    return None
