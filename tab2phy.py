#!/usr/bin/env python
import sys

import auto_phylip

if __name__ == '__main__':
    for file in sys.argv[1:]:
        auto_phylip.tab2phy(file)
