JumboDB tool for de Bruijn graph construction
==============

### Version: 1.0.1

This is jumboDB tool for fast de Bruijn graph construction from long sequences (reads or genomes) with very low error rate.
JumboDB is not a genome assembler by itself but rather a subroutine that translates a set of reads into compressed de Bruijn graph.
JumboDB a part of LJA (La Jolla Assembler) genome assembler under development that is designed for HiFi read assembly.
Unique feature of jumboDB is that it can construct de Bruijn graphs for any value of k.
Moreover increasing k does not lead to significant increase in time and space requirements.
Primarily time and space performance depend on the size of the resulting graph (total length of all edges in nucleotides).

JumboDB uses a combination of multiple known techniques such as bloom filters, sparse de Bruijn graphs and rolling hashs.
JumboDB constructs de Bruijn graph of human genome from HiFi dataset with coverage 25 within 4 hours and using 32 threads on any value of k from 250 to 7000.
For error-corrected reads only 30 minutes is enough.

Credits
-------

This tool is developed by Anton Bankevich in [Pavel Pevzner's lab at UCSD](http://cseweb.ucsd.edu/~ppevzner/).
Pease cite [this paper](https://www.biorxiv.org/content/10.1101/2020.12.10.420448v1) if you use jumboDB in your work.

