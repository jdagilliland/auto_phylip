#!/usr/bin/env python
import auto_phylip

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
        description=("Run PHYLIP program to build trees from sequences in a" +
            " PHYLIP formatted file."),
        )
    parser.add_argument('files', nargs='+')
    argspace = parser.parse_args()
    for fname in argspace.files:
        auto_phylip.run_phylip(fname)

