'''
auto_phylip contains functions to handle Phylip programs.
'''
import csv

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
    return None
