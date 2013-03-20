"""The **Wedring** pipeline module.
The **Wedring** pipeline gathers some bioinformatic softwares to achieve
differentially expressed genes from RNA-Seq experiments.

The softwares used are:
* `Bowtie <http://bowtie-bio.sourceforge.net/index.shtml>`_
* `TopHat <http://tophat.cbcb.umd.edu/>`_
* `SAMtools <http://samtools.sourceforge.net/>`_
* `BEDTools <http://code.google.com/p/bedtools/>`_
* `DESeq. <http://www-huber.embl.de/users/anders/DESeq/>`_

The pipeline has three stages:
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

.. moduleauthor:: Rafael Covre <covrebio06@gmail.com>

"""

__version__ = "0.3.2"

from wedring.gffutils import *
from wedring.indexstage import *
from wedring.manager import *
from wedring.mapstage import *
from wedring.sysutils import *
from wedring.tableutils import *
from wedring.wedrerror import *
