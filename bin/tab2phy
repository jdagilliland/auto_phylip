#!/usr/bin/env python
import auto_phylip

if __name__ == '__main__':
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
    auto_phylip.tab2phy(
        argspace.files,
        outfile=argspace.phyfname,
        match=argspace.match,
        column=argspace.column,
        )
