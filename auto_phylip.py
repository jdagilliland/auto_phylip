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
    n_clone = len(lst_seqpair)
    # WARNING here I assume that all sequences in the tabfile have the same
    # length
    len_clone = len(lst_seqpair[0][1])
    if outfile == None:
        outfile = tabfile.rpartition('.')[0] + '.phy'
    with open(outfile,'wb') as f:
        f.write(str(n_clone).ljust(10,' ') + str(len_clone) + '\n')
        for seqpair in lst_seqpair:
            f.write(seqpair[0].ljust(10,' ') + seqpair[1] + '\n')
    return None
