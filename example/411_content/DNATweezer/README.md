DNATweezer
==========

DNATweezer is a collection of wrapper scripts for various modules in [BioPerl](http://BioPerl.org).

Tools
-----
The tools in DNATweezer make it easy to access the functionality of some of the BioPerl modules from the command line, without the need to write your own scripts for small tasks. (Even some larger tasks!)

For example, you could remove a sequence with id "E. coli" in a FASTA file, simply by doing:
    seq-manipulations.pl -d 'id:E. coli' input_file > new_file

The tools included are:

- **`bioaln`**: A wrapper for Bio::SimpleAlign and Bio::AlignIO.
- **`bioseq`**: A wrapper for Bio::Seq and Bio::SeqIO.
- **`biopop`**: A wrapper for Bio::AlignDNAStatistics, Bio::PopGen::Utilities and other BioPerl population genetics modules.
- **`biotree`**: A wrapper for Bio::Tree::Node, Bio::Tree::Tree and Bio::TreeIO.

A few other minor tools are included as well.

**NOTE** Only `bioaln` and `bioseq` are considered production quality at the moment. The other tools are still under heavy development.

Support
-------
As with any software, these tools can contain bugs and are always going to be works in progress. Thus, feedback such as bug reports and feature requests are welcomed.
