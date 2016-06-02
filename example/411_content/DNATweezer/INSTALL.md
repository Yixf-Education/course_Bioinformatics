Instructions
==========

All executable files are found in the top level of this directory.

Currently, the two most mature tools are "bioseq" and "bioaln".


Installation
-----
Instructions apply to Linux or UNIX-based systems (Mac, BSD, etc).
Suggestions for instructions for Windows are welcome.

**Method 1.**

Simply decompress anywhere and add the directory containing all files to your `PATH`

**Method 2.**

1. Copy the files "`bioseq`" and "`bioaln`" to a directory in your path. (e.g., `/usr/local/bin`).

2. Copy the file "`lib/SeqManipulations.pm`" to a directory in Perl's `@INC` (e.g., `/usr/local/lib`).

Run the tools by typing in their names, e.g., `bioaln` or `bioseq`.

You can get help using them at any time by simply supplying the help or man options:
    bioaln --help # Print help. Also can use bioaln -h
    bioseq --man # Print the manual page
