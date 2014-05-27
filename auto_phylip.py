'''
auto_phylip contains functions to handle Phylip programs.
'''
import csv
import os
import subprocess as sub
import re

phy_exec_default = ['phylip', 'dnapars']
boot_exec_default = ['phylip', 'seqboot']
cons_exec = ['phylip', 'consense']

n_bootstrap_default = 1000

def tab2phy(lst_tabfile, germline=None, outfile=None, **kwarg):
    """
    Generate a PHYLIP formatted *.phy file from a list of tabfiles.

    A column and a regex may be supplied which will be used to filter
    the results from the tabfile.

    Parameters
    ----------
    lst_tabfile : list
        A list of tabfiles from which to draw entries to populate the
        *.phy file.
    column : str
        The name of the column whose value must match `match` in order
        for a given row to be included in the output *.phy file.
            (default: 'CLONE')
    match : str
        A string which will be converted to a regex to match to entries
        in the tabfiles. (default: None; all rows will be included)
    flags : int
        This is a sum of flags to be passed to the regex compiler.
        (default: 0; default regex settings)
    """
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
    """
    Convert a tabfile entry into a pair specifying the sequence name to
    be used in the *.phy file and the sequence itself.
    """
    return (entry['SEQUENCE_ID'][-9:], entry['SEQUENCE'])

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
    with open(outfile, 'wb') as f:
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
    with open(tabfile, 'rb') as f:
        reader = csv.DictReader(f, delimiter='\t')
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
    # elif type(match) == _sre.SRE_Pattern:
    #     # check if the provided match is a compiled regex, if so, use as is
    #     reg_match = match
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
    """
    Generate a string from a name and a sequence which can be written as
    a line in a *.phy file.
    """
    return str(name).ljust(10, ' ') + str(sequence) + '\n'

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
    # if bootstrap is provided use, otherwise None
    bootstrap = kwarg.pop('bootstrap', None)
    if bootstrap:
        # save the original input *.phy name
        phy_in_orig = phy_in
        # for most of following operations, use bootstrapped *.phy file
        phy_in = run_seqboot(phy_in, bootstrap, **kwarg)
    basename = phy_in.rpartition('.')[0]
    # remove old files
    if os.path.lexists('infile'):
        os.remove('infile')
    if os.path.lexists('outfile'):
        os.remove('outfile')
    if os.path.lexists('outtree'):
        os.remove('outtree')
    # write a command file with the specified options
    lst_phy_opts = _get_phy_opts(phy_in, bootstrap=bootstrap, **kwarg)
    cmdfname = write_cmdfile(lst_phy_opts, trailing_nl=False)
    # open phylip process
    p = sub.Popen(phy_exec, stdin=sub.PIPE)
    # send command file contents as input, wait for output
    p.communicate(open(cmdfname, 'r').read())
    # rename output files
    try:
        os.rename('outfile', basename + '.out')
        os.rename('outtree', basename + '.tree')
    except:
        print('The expected output was not generated. Phylip may have failed')
        raise
    try:
        os.remove(cmdfname)
    except:
        raise
    if bootstrap:
        # run consense only if bootstrapping was performed
        run_consense(basename + '.tree', **kwarg)
        print('Consensus tree is: {name}'.format(name=basename + '.tree'))
        print('Using original phy file: {name}'.format(name=phy_in_orig))
        cleanup_consense(basename + '.tree', phy_in_orig)
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
    p.communicate(open(cmdfname, 'r').read())
    basename = fname.rpartition('.')[0]
    bootname = basename + '.boot.phy'
    try:
        os.rename('outfile', bootname)
    except:
        print('The expected output was not generated. Phylip may have failed')
        raise
    try:
        os.remove(cmdfname)
    except:
        raise
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
    p.communicate(open(cmdfname, 'r').read())
    basename = fname.rpartition('.')[0]
    try:
        os.rename('outfile', basename + '.cons.out')
        os.rename('outtree', basename + '.cons.tree')
    except:
        print('The expected output was not generated. Phylip may have failed')
        raise
    try:
        os.remove(cmdfname)
    except:
        raise
    return None

def cleanup_consense(fname_consensus, fname_orig, **kwarg):
    """
    Use a consensus tree, and an original (non-bootstrapped) *.phy input
    file to generate a cleaned up consensus tree that estimates the true
    tree.
    """
    # if phy_exec is provided, or None, use default
    phy_exec = kwarg.pop('phy_exec', None)
    if phy_exec == None:
        phy_exec = phy_exec_default
    basename = fname_orig.rpartition('.')[0]
    # remove old files
    if os.path.lexists('infile'):
        os.remove('infile')
    if os.path.lexists('intree'):
        os.remove('intree')
    if os.path.lexists('outfile'):
        os.remove('outfile')
    if os.path.lexists('outtree'):
        os.remove('outtree')
    # setup options for consensus cleanup
    lst_phy_opts = _get_phy_opts(fname_orig,
            fname_tree=fname_consensus,
            search=False,
            )
    cmdfname = write_cmdfile(lst_phy_opts, trailing_nl=False)
    # open phylip process
    p = sub.Popen(phy_exec, stdin=sub.PIPE)
    # send command file contents as input, wait for output
    p.communicate(open(cmdfname, 'r').read())
    # rename output files
    try:
        os.rename('outfile', basename + '.out')
        os.rename('outtree', basename + '.tree')
    except:
        print('The expected output was not generated. Phylip may have failed')
        raise
    try:
        os.remove(cmdfname)
    except:
        raise
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
    with open(cmdfname, 'w') as f:
        for arg in opts[:-1]:
            f.write(arg + '\n')
        f.write(opts[-1])
    return cmdfname

def _get_phy_opts(fname, **kwarg):
    """
    Based on the filename and optional arguments, get a list of options
    which can be written on lines of a command file to be fed into
    PHYLIP phylogeny programs.
    """
    # If seed is not specified, use 9
    seed = kwarg.pop('seed', 9) # must be 4n+1 for some reason!?
    # If jumble is not specified, do not jumble
    jumble = kwarg.pop('jumble', 1)
    bootstrap = kwarg.pop('bootstrap', False)
    # if bootstrap:
    #     seed = True
    search = kwarg.pop('search', True)
    fname_tree = kwarg.pop('fname_tree', None)
    opts = list()
    # first the filename
    opts.append(fname)
    # set search option
    #-----------------
    if search:
        opts.append('S')#|
        opts.append('Y')#|
    else:
        opts.append('U')
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
    if not search:
        # Phylogeny programs prompt for the filename of the user tree if
        # told not to search for the best tree.
        opts.append(fname_tree)
        opts.append(str(seed))
    return opts

def _get_seqboot_opts(fname, n_bootstrap, weights=False, seed=9, **kwarg):
    """
    Based on the filename, bootstrap number, and optional arguments,
    get a list of options which can be written on lines of a command
    file to be fed into PHYLIP seqboot.
    """
    lst_seqboot_opts = []
    lst_seqboot_opts.append(fname)
    lst_seqboot_opts.append('R')
    lst_seqboot_opts.append(str(n_bootstrap))
    # the following to be properly implemented later
    if weights:
        lst_seqboot_opts.append('W')
    # finish putting in options, and confirm
    lst_seqboot_opts.append('Y')
    lst_seqboot_opts.append(str(seed))
    return lst_seqboot_opts

def _get_consense_opts(fname, **kwarg):
    """
    Based on the filename, and optional arguments,
    get a list of options which can be written on lines of a command
    file to be fed into PHYLIP consense.
    """
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
    """
    The main runner script for the command `run_phylip`, including the
    argparse parser.
    """
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
            NOTE: Use of this option is not recommended by the author,
            use at your own risk.
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
    parser.add_argument('-s', '--seed',
            dest='seed',
            default=9,
            type=int,
            help="""
            Random seed to use for some PHYLIP programs.
            """,
            )
    parser.add_argument('-j', '--jumble',
            dest='jumble',
            default=1,
            type=int,
            help="""
            The number of times to jumble the order of the input
            sequences.
            The more you jumble, the longer run times you will have
            (about linear), but it does seem to prevent sub-optimal
            trees coming out.
            """
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
                seed=argspace.seed,
                jumble=argspace.jumble,
                )

def _run_seqboot_main():
    """
    The main runner script for the command `run_seqboot`, including the
    argparse parser.
    """
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
    parser.add_argument('-s', '--seed',
            dest='seed',
            default=9,
            type=int,
            help="""
            Random seed to use for seqboot.
            """,
            )
    parser.add_argument('files', nargs='+',
            help="""
            These are the PHY files from which to bootstrap expanded
            datasets.
            They are treated independently.
            """,
            )
    argspace = parser.parse_args()
    for fname in argspace.files:
        run_seqboot(fname, argspace.bootstrap,
                seed=argspace.seed,
                )
    return None

def _run_consense_main():
    """
    The main runner script for the command `run_consense`, including the
    argparse parser.
    """
    import argparse
    parser = argparse.ArgumentParser(
            description='''Build a consensus tree from bootstrapped trees''',
            )
    parser.add_argument('files', nargs='+',
            help="""
            These are the file[s] from which to build consensus trees.
            Each file should be a self-contained set of trees delimited
            by semicolons (;).
            The provided files are treated independently.
            """,
            )
    argspace = parser.parse_args()
    for fname in argspace.files:
        run_consense(
                fname,
                )
    return None

def _cleanup_consense_main():
    """
    The main runner script for the command `cleanup_consense`, including the
    argparse parser.
    """
    import argparse
    parser = argparse.ArgumentParser(
            description='''Correct the branch lengths on a consensus tree''',
            )
    parser.add_argument('-p', '--phyfile',
            dest='phyfile',
            help="""
            This PHY file will be used to determine the length of the
            sequences in the tree, and that way correct the branch
            lengths.
            """,
            )
    parser.add_argument('consensus', nargs='+',
            help="""
            This is the consensus tree to clean up.
            Clean up is done in place, so if you want to preserve the
            original state of the tree prior to cleanup, you should copy
            it.
            You can do this to multiple trees at once, as long as they
            use the same PHY file, i.e. have the same sequence length.
            """,
            )
    argspace = parser.parse_args()
    for fname in argspace.consensus:
        cleanup_consense(fname, argspace.phyfile,)
    return None

def _tab2phy_main():
    """
    The main runner script for the command `tab2phy`, including the
    argparse parser.
    """
    import argparse
    parser = argparse.ArgumentParser(
        description='Build a phylip formatted file from a tabfile.',
        )
    parser.add_argument('files', nargs='+',
        help="""
        These are the TAB file[s] that you provide which you use to
        build the PHY file.
        Order is important, since some downstream utilities which use
        the PHY file will treat the first sequence therein as a root
        sequence.
        """,
        )
    parser.add_argument('-m', '--match',
        dest='match',
        default=None,
        help="""
        Provide a regex to which to match the specified column of the
        TAB file.
        This should be a python style regex, since it goes directly into
        the python re module.
        """,
        )
    parser.add_argument('-f', '--field',
        dest='column',
        default='CLONE',
        help="""
        Specify the column whose contents to match according to the
        given regex.
        """,
        )
    parser.add_argument('-c', '--combine',
        dest='phyfname',
        nargs='?',
        const='file.phy',
        default=None,
        help="""
        Specify a filename to combine multiple TAB files into a single
        PHY file.
        Generally a better time to combine sequence data like this is
        when generating the tab files in the first place.
        """,
        )
    argspace = parser.parse_args()
    tab2phy(
        argspace.files,
        outfile=argspace.phyfname,
        match=argspace.match,
        column=argspace.column,
        )
    return None
