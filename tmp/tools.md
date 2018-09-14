# Tools

## [Converting Genome Coordinates From One Genome Version To Another](https://www.biostars.org/p/65558/)

Some recent posts reminded me that it might be useful for us to review the options for converting between genome coordinate systems.

This comes up in several contexts. Probably the most common is that you have some coordinates for a particular version of a reference genome and you want to determine the corresponding coordinates on a different version of the reference genome for that species. For example, you have a bed file with exon coordinates for human build GRC37 (hg19) and wish to update to GRCh38. By the way, for a nice summary of genome versions and their release names refer to the [Assembly Releases and Versions FAQ](http://genome.ucsc.edu/FAQ/FAQreleases.html)

Or perhaps you have coordinates of a gene and wish to determine the corresponding coordinates in another species. For example, you have coordinates of a gene in human GRCh38 and wish to determine corresponding coordinates in mouse mm10.

Finally you may wish to convert coordinates between coordinate systems within a single assembly. For example, you have the coordinates of a series of exons and you want to determine the position of these exons with respect to the transcript, gene, contig, or entire chromosome.

There are now several well known tools that can help you with these kinds of tasks:

1. **UCSC liftOver**. This tool is available through a simple [web interface](http://genome.ucsc.edu/cgi-bin/hgLiftOver) or it can be downloaded as a [standalone executable](http://hgdownload.cse.ucsc.edu/admin/exe/). To use the executable you will also need to download the appropriate [chain file](http://hgdownload.cse.ucsc.edu/downloads.html#liftover). Each chain file describes conversions between a pair of genome assemblies. Liftover can be used through [Galaxy](https://usegalaxy.org/) as well. There is a python implementation of liftover called [pyliftover](https://pypi.python.org/pypi/pyliftover) that does conversion of point coordinates only.

2. **NCBI Remap**. This tool is conceptually similar to liftOver in that in manages conversions between a pair of genome assemblies but it uses different methods to achieve these mappings. It is also available through a simple [web interface](http://www.ncbi.nlm.nih.gov/genome/tools/remap) or you can use the [API for NCBI Remap](http://www.ncbi.nlm.nih.gov/genome/tools/remap/docs/api).

3. **The Ensembl API**. The final example I described above (converting between coordinate systems within a single genome assembly) can be accomplished with the [Ensembl core API](http://ensembl.org/info/docs/api/core/index.html#api). Many examples are provided within the [installation](http://ensembl.org/info/docs/api/api_installation.html), [overview](http://ensembl.org/info/docs/api/core/core_API_diagram.html), [tutorial](http://ensembl.org/info/docs/api/core/core_tutorial.html) and [documentation](http://ensembl.org/info/docs/Doxygen/core-api/index.html) sections of the Ensembl API project. In particular, refer to these sections of the [tutorial](http://ensembl.org/info/docs/api/core/core_tutorial.html): 'Coordinates', 'Coordinate systems', 'Transform', and 'Transfer'.

4. **Assembly Converter**. Ensembl also offers their own simple web interface for coordinate conversions called the [Assembly Converter](https://www.ensembl.org/Homo_sapiens/Tools/AssemblyConverter?db=core).

5. **Bioconductor rtracklayer** package. For [R](http://cran.us.r-project.org/index.html) users, [Bioconductor](http://www.bioconductor.org/) has an implementation of UCSC liftOver in the [rtracklayer package](http://bioconductor.org/packages/release/bioc/html/rtracklayer.html). To see documentation on how to use it, open an R session and run the following commands.
```
source("http://bioconductor.org/biocLite.R")
biocLite("rtracklayer")
library(rtracklayer)
?liftOver
```
6. **CrossMap**. A standalone [open source program](http://sourceforge.net/projects/crossmap/files/CrossMap-0.1.3.tar.gz/download) for convenient conversion of genome coordinates (or annotation files) between different assemblies. It supports most commonly used file formats including [SAM](http://samtools.sourceforge.net/SAM1.pdf)/BAM, Wiggle/BigWig, BED, GFF/GTF, VCF. CrossMap is designed to liftover genome coordinates between assemblies. It’s not a program for aligning sequences to reference genome. Not recommended for converting genome coordinates between species.

7. `Flo.` A liftover pipeline for different reference genome builds of the same species. It describes the process as follows: "align the new assembly with the old one, process the alignment data to define how a coordinate or coordinate range on the old assembly should be transformed to the new assembly, transform the coordinates."

8. [Picard Liftover VCF](https://broadinstitute.github.io/picard/command-line-overview.html#LiftoverVcf). Lifts over a VCF file from one reference build to another. This tool adjusts the coordinates of variants within a VCF file to match a new reference. The tool is based on the UCSC liftOver and uses a UCSC chain file to guide its operation.

