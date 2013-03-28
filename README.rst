==========================================================================================
**Wedring** - A pipeline for differentialy expressed genes analysis in RNA-Seq experiments
==========================================================================================

The **Wedring** pipeline gathers some bioinformatic softwares to achieve
differentially expressed genes from RNA-Seq experiments. It is implemented
as a Python package, named *wedring*, and also has the command line interface
*wedr*.

The softwares used are:
* `Bowtie <http://bowtie-bio.sourceforge.net/index.shtml>`_
* `TopHat <http://tophat.cbcb.umd.edu/>`_
* `SAMtools <http://samtools.sourceforge.net/>`_
* `BEDTools <http://code.google.com/p/bedtools/>`_
* `DESeq <http://www-huber.embl.de/users/anders/DESeq/>`_

The pipeline has four stages:
Pre-processing stage:
    The annotation file (in GFF format) is validated, according to its field numbers
    and new-line character.
Indexing stage:
    The Burrows-Wheeler (BW) index is created from a reference genome using the
    *bowtie-build* program. This stage is optional since a prebuilt index may
    be provided to the pipeline.
Mapping stage:
    RNA-Seq reads from different experimental conditions are mapped against the
    BW index. The mapping program may be *Bowtie* or *TopHat*, according to the
    user needing. *Bowtie* is recommended to analyze prokaryotic data.
    Meanwhile, *TopHat* is recommended to eukaryotic data due to the splicing
    junctions routines in its pipeline. Afterwards, informations from the
    mappings are processed using *SAMtools* and *BEDTools* to generate a
    counting table of each genomic features in each experimental condition.
Differential expression stage:
    The counting table is used as the input to the DESeq pipeline. The result
    will be a table containing statistics of each genomic feature which
    includes *p* and *fold change* values.