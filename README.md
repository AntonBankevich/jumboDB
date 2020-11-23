JumboDB tool for de Bruijn graph construction
==============

### Version: 1.0

This is jumboDB tool for fast de Bruijn graph construction from long sequences (reads or genomes) with very low error rate.
JumboDB is a part of LJA (La Jolla Assembler) genome assembler for HiFi reads.
Unique feature of jumboDB is that it can construct de Bruijn graphs for any value of k.
Moreover increasing k does not lead to significant increase in time and space requirements.
Primarily time and space performance depend on the size of the resulting graph.

Note that jumboDB discards all reads that are shorter than k + w (see [manual](manual.md)).

JumboDB uses a combination of multiple known techniques such as bloom filters, sparse de Bruijn graphs and rolling hashs.
JumboDB constructs de Bruijn graph of human genome from HiFi dataset with coverage 25 within 4 hours and using 32 threads on any value of k from 250 to 7000.
For error-corrected reads only 30 minuts is enough.

For installation and running instructions please refer to [manual](manual.md)

License
-------

This tool is distributed under a BSD license. See the [LICENSE file](LICENSE) for details.


Credits
-------

This tool is developed by Anton Bankevich in [Pavel Pevzner's lab at UCSD](http://cseweb.ucsd.edu/~ppevzner/)
