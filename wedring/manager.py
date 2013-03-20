"""The **Wedring** pipeline manager module.
This module provides the :class:Wedring class which is the manager of the
pipeline. It executes all the pipeline steps.

.. moduleauthor:: Rafael Covre <covrebio06@gmail.com>

"""

from getopt import getopt, GetoptError
from itertools import count
from multiprocessing import Pool
from os import environ, path
from types import MethodType
try:
    from ConfigParser import RawConfigParser
except ImportError:
    from configparser import RawConfigParser
try:
    from copy_reg import pickle
except ImportError:
    from copyreg import pickle

from wedring.gffutils import GffFormatError, write_validated_gff
from wedring.indexstage import IndexBuilder
from wedring.mapstage import WedringMast
from wedring.sysutils import (BioSoft, wedr_check_path, wedr_prepare_directory,
                              wedr_clean, wedr_report, wedr_which)
from wedring.tableutils import (write_genomic_features_to_file, prepare_table_header,
                                write_count_table_to_file)
from wedring.wedrerror import WedringError

__all__ = ["Wedring"]

TTL_PIPELINE = 0
JUST_INDEX = 1
JUST_MAP = 2
JUST_TABLE = 3
JUST_DE = 4

def _method_pickle(method):
    return _method_unpickle, (method.im_func.__name__,
                              method.im_self,
                              method.im_class)

def _method_unpickle(fx_name, obj, cls):
    for c in cls.mro():
        try:
            fx = c.__dict__[fx_name]
        except KeyError:
            pass
        else:
            break
    return fx.__get__(obj, cls)

pickle(MethodType, _method_pickle, _method_unpickle)

class Wedring(object):
    """This class is the the manager of the **Wedring** pipeline. It manages all
    the parameters and all the steps of the pipeline.

    """

    def __init__(self, arguments=None):
        """:class:Wedring constructor.

        :param arguments: List of arguments, usually command line ones.
        :type arguments: list

        """
        self._wedr_list = []
        self._num_threads = 1
        self._qt = False             # Quietly?
        self._o = "./wedr_out"       # Output directory.
        self._id = "./wedr_index"    # Bowtie's BW index output directory.
        self._m = "bowtie"           # Mapper.
        self._r = None               # Reference sequence.
        self._i = None               # Bowtie's EBWT index.
        self._l = [None]             # Single-end library.
        self._1 = [None]             # Pair-end library mate 1.
        self._2 = [None]             # Pair-end library mate 2.
        self._q = [None]             # Single-end library quality.
        self._q1 = [None]            # Pair-end library mate 1 quality.
        self._q2 = [None]            # Pair-end library mate 2 quality.
        self._cf = [None]            # Coverage files.
        self._a = None               # Annotation file.
        self._c = None               # Configuration file.
        self._ml = [None]            # Mapping label.
        self._il = None              # Index label.
        self._cnd = None             # Experiment conditions.
        self._cnt_table = None       # Name of the counts table.
        self._ib = False             # Build Bowtie's BW index?
        self._wb = TTL_PIPELINE      # Wedring's barrier.
        self._indexbldr = None       # Wedring's IndexBuilder instance.
        self._ld = None              # Directory where log files are saved.
        self.wedring_parse_args(arguments)

    def wedring_parse_args(self, arguments):
        """This method parses the :class:Wedring arguments and make some
        adjustments.

        :param arguments: List of arguments, usually command line ones.
        :type arguments: list
        :raises: :class:WedringError

        """
        try:
            if len(arguments) < 1:
                raise WedringError(133, "Insufficient number of arguments.")
        except TypeError:
            raise WedringError(133, "Arguments were not provided.")
        bin_path = None # path to add to system path

        group_all = False
        if "--group-all" in arguments:
            group_all = True
            arguments.remove("--group-all")

        # Lambda function to adjust parameters according to the --group-all
        # argument value.
        adjust_param = lambda params: ([param for param in params.split(',')
                                        if param != '']
                                       if not group_all else [param])

        # Adjusting the value of the Wedring barrier to define wich pipeline
        # stage will be executed
        if "--just-indexbuild" in arguments:
            self._wb = JUST_INDEX
            arguments.remove("--just-indexbuild")
        if "--just-map" in arguments:
            self._wb = JUST_MAP
            arguments.remove("--just-map")
        if "--just-counttable" in arguments:
            self._wb = JUST_TABLE
            arguments.remove("--just-counttable")
        if "--just-de" in arguments:
            self._wb = JUST_DE
            arguments.remove("--just-de")

        try:
            opts = getopt(arguments, "n:o:x:m:r:i:l:1:2:q:a:g:t:c:t:d:p:",
                          ["num-threads=",
                           "out-dir=",
                           "index-dir=",
                           "quiet",
                           "mapper=",
                           "ref-sequence=",
                           "bw-index=",
                           "lib-file=",
                           "pair-mate-1=",
                           "pair-mate-2=",
                           "quals=",
                           "q1=",
                           "q2=",
                           "annot-file=",
                           "coverage-files=",
                           "config-file=",
                           "count-table=",
                           "map-label=",
                           "index-label=",
                           "conditions=",
                           "path="])
            if opts[1] != []:
                raise WedringError(132, "Argument list not supported: %s." %
                                   " ".join(opts[1]))
            for opt, val in opts[0]:
                if opt in ("-n", "--num-threads"):
                    self._num_threads = int(val)
                elif opt in ("--quiet"):
                    self._qt = True
                elif opt in ("-m", "--mapper"):
                    self._m = val
                elif opt in ("-r", "--ref-sequence"):
                    self._r = val
                elif opt in ("-i", "--bw-index"):
                    self._i = val
                elif opt in ("-a", "--annot-file"):
                    self._a = val
                elif opt in ("-c", "--config-file"):
                    self._c = val
                elif opt in ("-t", "--count-table"):
                    self._cnt_table = val
                elif opt == "--index-label":
                    self._il = val
                elif opt in ("-d", "--conditions"):
                    self._cnd = val
                elif opt in ("-x", "--index-dir"):
                    self._id = val
                elif opt in ("-o", "--out-dir"):
                    self._o = val
                elif opt in ("-l", "--lib-file"):
                    self._l = adjust_param(val)
                elif opt in ("-1", "--pair-mate-1"):
                    self._1 = adjust_param(val)
                elif opt in ("-2", "--pair-mate-2"):
                    self._2 = adjust_param(val)
                elif opt in ("-q", "--quals"):
                    self._q = adjust_param(val)
                elif opt == "--q1":
                    self._q1 = adjust_param(val)
                elif opt == "--q2":
                    self._q2 = adjust_param(val)
                elif opt in ("-g", "--coverage-files"):
                    if ',' not in val:
                        raise WedringError(134, "You must provide a list of coverage files.")
                    else:
                        self._cf = [cf for cf in val.split(',') if cf != '']
                elif opt == "--map-label":
                    self._ml = adjust_param(val)
                elif opt in ("-p", "--path"):
                    bin_path = val
        except GetoptError as err:
                raise WedringError(136, "%s." % str(err).capitalize())
        if bin_path is not None:
            environ["PATH"] = path.pathsep.join((environ["PATH"], path))
        if self._wb != JUST_INDEX:
            self._ld = path.join(self._o, "log")
            if self._a:
                if not self._qt:
                    wedr_report("Validating GFF file: \'%s\'." % self._a)
                try:
                    gff_out = write_validated_gff(self._a)
                except GffFormatError as gffferr:
                    raise WedringError(134, "[%s] %s" %
                                        (type(gffferr).__name__, gffferr))
                except EnvironmentError as env_err:
                    raise WedringError(env_err.errno, "[%s (%d)] %s%s%s." %
                                (type(env_err).__name__,
                                 env_err.errno, 
                                 env_err.strerror,
                                 ': ' if env_err.filename else '',
                                 env_err.filename if env_err.filename else ''))
                if isinstance(gff_out, str):
                    if not self._qt:
                        wedr_report("Now using validated GFF file: \'%s\'." %
                                     gff_out)
                    self._a = gff_out
                    if not self._qt:
                        wedr_report("Gff validation - Done!")
                else:
                    if not self._qt:
                        wedr_report("Gff validation - Done!")
            if self._wb in (TTL_PIPELINE, JUST_TABLE, JUST_DE):
                if self._cnd is None:
                    raise WedringError(135, "You must set the experimental conditions.")
            if self._wb in (TTL_PIPELINE, JUST_MAP):
                if self._l != [None] and self._q == [None]:
                    self._q += [None] * (len(self._l) - 1)
                if self._l != [None] and self._ml == [None]:
                    self._ml += [None] * (len(self._l) - 1)
                if (self._1 != [None] and self._2 != [None] and
                    self._q1 == [None] and self._q1 == [None]):
                    self._q1 += [None] * (len(self._1) - 1)
                    self._q2 += [None] * (len(self._2) - 1)
                if self._1 != [None] and self._2 != [None] and self._ml == [None]:
                    self._ml += [None] * (len(self._1) - 1)
            elif self._wb == JUST_DE and self._cnt_table is None:
                raise WedringError(135, "You must provide the counting table.")

    def wedring_indexbuilder_parse_params(self):
        """This method sets **Wedring**'s index builder specific parameters."""
        if self._i is None:
            self._ib = True
            self._indexbldr = IndexBuilder(out_dir=self._id,
                                           ref_seq=self._r,
                                           index_label=self._il,
                                           cfg_file=self._c,
                                           quiet=self._qt)

    def wedring_mapping_parse_params(self):
        """This method sets **Wedring**'s mapping stage specific parameters.

        :raises: :class:WedringError

        """
        wedr_append = self._wedr_list.append
        if self._cnd is None:
            _next = count(1).next
        else:
            _next = prepare_table_header(self._cnd).next
        if self._wb != JUST_TABLE:
            if self._wb != JUST_INDEX:
                if self._l != [None]:
                    if self._q != [None]:
                        for idx, lib in enumerate(self._l):
                            try:
                                wedr = WedringMast(mapper=self._m,
                                                   out_dir=self._o,
                                                   index=self._i,
                                                   lib_file=lib,
                                                   qual_file=self._q[idx],
                                                   annot_file=self._a,
                                                   cfg_file=self._c,
                                                   map_label=self._ml[idx],
                                                   barrier=self._wb,
                                                   quiet=self._qt,
                                                   id_=_next())
                            except IndexError:
                                raise WedringError(140, "Unbalanced options -l/--lib-file, --quals, --map-label or -o/--out-dir.")
                            wedr_append(wedr)
                    else:
                        for idx, lib in enumerate(self._l):
                            try:
                                wedr = WedringMast(mapper=self._m,
                                                   out_dir=self._o,
                                                   index=self._i,
                                                   lib_file=lib,
                                                   annot_file=self._a,
                                                   cfg_file=self._c,
                                                   map_label=self._ml[idx],
                                                   barrier=self._wb,
                                                   quiet=self._qt,
                                                   id_=_next())
                            except IndexError:
                                raise WedringError(140, "Unbalanced options -l/--lib-file, --map-label or -o/--out-dir.")
                            wedr_append(wedr)
                elif self._1 != [None]:
                    if self._q1 != [None]:
                        for idx, lib_1 in enumerate(self._1):
                            try:
                                wedr = WedringMast(mapper=self._m,
                                                   out_dir=self._o,
                                                   index=self._i,
                                                   lib_mate_1=lib_1,
                                                   lib_mate_2=self._2[idx],
                                                   q_mate_1=self._q1[idx],
                                                   q_mate_2=self._q2[idx],
                                                   annot_file=self._a,
                                                   cfg_file=self._c,
                                                   map_label=self._ml[idx],
                                                   barrier=self._wb,
                                                   quiet=self._qt,
                                                   id_=_next())
                            except IndexError:
                                raise WedringError(140, "Unbalanced options -1/--pair-mate-1, -2/--pair-mate-2, --q1, --q2, --map-label or -o/--out-dir.")
                            wedr_append(wedr)
                    else:
                        for idx, lib_1 in enumerate(self._1):
                            try:
                                wedr = WedringMast(mapper=self._m,
                                                   out_dir=self._o,
                                                   index=self._i,
                                                   lib_mate_1=lib_1,
                                                   lib_mate_2=self._2[idx],
                                                   annot_file=self._a,
                                                   cfg_file=self._c,
                                                   map_label=self._ml[idx],
                                                   barrier=self._wb,
                                                   quiet=self._qt,
                                                   id_=_next())
                            except IndexError:
                                raise WedringError(140, "Unbalanced options -1/--pair-mate-1, -2/--pair-mate-2, --map-label or -o/--out-dir.")
                            wedr_append(wedr)
        else:
            for cf in self._cf:
                wedr_append(WedringMast(annot_file=self._a,
                                        cov_file=cf,
                                        barrier=self._wb,
                                        id_=_next()))

    def exec_deseq(self):
        """Execute the **Wedring**'s differential expression stage.

        :raises: :class:WedringError

        """
        if not self._qt:
            wedr_report("Calculating differential expression with DESeq.")
        outfile = path.join(self._o, "diffexpr.txt")
        errfile = path.join(self._ld, "diffexpr.log")
        de_cmd_line = "Rscript --vanilla %s %s %s %s %s" % (wedr_which("diffExprStage.R"),
                                                            self.wedring_diffexpr_parse_options(),
                                                            self._cnt_table,
                                                            self._cnd,
                                                            outfile)
        de = BioSoft(de_cmd_line, errfile=errfile)
        if not self._qt:
            wedr_report("Command line:\n    %s" % de.command)
        de.run()
        if 0 != de.return_code != -1:
            raise WedringError(141, "%s exitted with status %d. See log file '%s' for more details." %
                               (de.program_name, de.return_code, de.errfile))
        wedr_clean(de.errfile)
        # TODO Add verification of the DESeq's output with wedr_check_path()
        #   \_ table (OK), graphics
        wedr_check_path(outfile)
        if not self._qt:
            wedr_report("DESeq - Done!")

    def wedring_diffexpr_parse_options(self):
        """Parser for the DE section of the configuration file.

        :returns: Additional options for the **Wedring**'s Differential
        Expression stage.

        """
        de_cmd = ""
        params = {"method": "pooled",
                  "sharing_mode": "maximum",
                  "fit_type": "parametric",
                  "ma_plot": "false",
                  "volcano_plot": "false",
                  "dispest_plot": "false",
                  "alpha": "0.05",
                  "img_size": "medium"}
        cf_parser = RawConfigParser()
        cf_parser.read(self._c)
        for param, val in cf_parser.items("DE"):
            params[param] = val
        for param, val in params.iteritems():
            if param == "method" and val != "pooled":
                de_cmd += " -method=%s" % val
            elif param == "sharing_mode" and val != "maximum":
                de_cmd += " -sharing.mode=%s" % val
            elif param == "fit_type" and val != "parametric":
                de_cmd += " -fit.type=%s" % val
            elif param == "ma_plot" and val != "false":
                de_cmd += " -ma.plot"
            elif param == "volcano_plot" and val != "false":
                de_cmd += " -volcano.plot"
            elif param == "dispest_plot" and val != "false":
                de_cmd += " -dispest.plot"
            elif param == "p_value" and val != "0.05":
                de_cmd += " -p.value=%s" % val
            elif param == "img_size" and val != "medium":
                de_cmd += " -img.size=%s" % val
        return de_cmd

    def run(self):
        """Execute all steps of the **Wedring** pipeline."""
        # This method will execute according to the self._wb value:
        # The value are set after the command line options --just-indexbuild,
        # --just-map, --just-counttable, --just-de, and the possible values are
        # defined in the globals TTL_PIPELINE, JUST_INDEX, JUST_MAP, JUST_TABLE
        # and JUST_DE, which mean:
        # TTL_PIPELINE -- execute all steps of the pipeline
        # JUST_INDEX -- just execute the indexing stage
        # JUST_MAP -- execute the indexing stage (if needed) and the mapping
        #             stage
        # JUST_TABLE -- just build the count table
        # JUST_DE -- just execute the differential expression stage
        if self._wb in (TTL_PIPELINE, JUST_INDEX, JUST_MAP):
            self.wedring_indexbuilder_parse_params()
            if self._indexbldr is not None:
                self._indexbldr.run()
                self._i = self._indexbldr.index
            if self._wb != JUST_INDEX:
                wedr_prepare_directory(self._o)
                wedr_prepare_directory(self._ld)
                self.wedring_mapping_parse_params()
                p = Pool(self._num_threads)
                self._wedr_list = p.map(WedringMast.exec_mapping_stage,
                                        self._wedr_list)
                if self._wb != JUST_MAP:
                    feats_file = path.join(self._o, "genomic_features.txt")
                    tbl_file = path.join(self._o, "count_table.txt")
                    if not self._qt:
                        wedr_report("Writing genomic features to disk.")
                    write_genomic_features_to_file(self._a, feats_file)
                    if wedr_check_path(feats_file):
                        if not self._qt:
                            wedr_report("Writing genomic features - Done!")
                    cov_f = [wedrmast.cov_file for wedrmast in self._wedr_list]
                    if not self._qt:
                        wedr_report("Writing count table to disk.")
                    write_count_table_to_file(self._a, cov_f, self._cnd,
                                               tbl_file)
                    if wedr_check_path(tbl_file):
                        self._cnt_table = tbl_file
                        if not self._qt:
                            wedr_report("Writing count table - Done!")
                    self.exec_deseq()
        elif self._wb == JUST_TABLE:
            self.wedring_mapping_parse_params()
            wedr_prepare_directory(self._o)
            feats_file = path.join(self._o, "genomic_features.txt")
            tbl_file = path.join(self._o, "count_table.txt")
            if not self._qt:
                wedr_report("Writing genomic features to disk.")
            write_genomic_features_to_file(self._a, feats_file)
            if wedr_check_path(feats_file):
                if not self._qt:
                    wedr_report("Writing genomic features - Done!")
            cov_f = [wedrmast.cov_file for wedrmast in self._wedr_list]
            if not self._qt:
                wedr_report("Writing count table to disk.")
            write_count_table_to_file(self._a, cov_f, self._cnd, tbl_file)
            if wedr_check_path(tbl_file):
                self._cnt_table = tbl_file
                if not self._qt:
                    wedr_report("Writing count table - Done!")
        elif self._wb == JUST_DE:
            wedr_prepare_directory(self._o)
            wedr_prepare_directory(self._ld)
            self.exec_deseq()
