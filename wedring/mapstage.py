"""The **Wedring** pipeline mapping stage module.
This module provides the :class:WedringMast class which comprises the mapping
stage of the **Wedring** pipeline.

.. moduleauthor:: Rafael Covre <covrebio06@gmail.com>

"""

from os import path as path
try:
    from ConfigParser import RawConfigParser
except ImportError:
    from configparser import RawConfigParser

from wedring.sysutils import (BioSoft, wedr_check_path, wedr_check_program,
                              wedr_clean, wedr_report)
from wedring.wedrerror import WedringError

__all__ = ["WedringMast"]

class WedringMast(object):
    """**Wedring** pipeline mapping stage."""

    def __init__(self, out_dir="./wedr_out", mapper="bowtie", index=None,
                 lib_file=None, lib_mate_1=None, lib_mate_2=None,
                 qual_file=None, q_mate_1=None, q_mate_2=None, annot_file=None,
                 cov_file=None, cfg_file=None, map_label=None, barrier=0,
                 quiet=False, id_=1):
        """:class:WedringMast constructor.

        :param out_dir: Which output directory? (default: ./wedr_out)
        :type out_dir: str -- Pathname to an existing (or not) directory.
        :param mapper: Which mapper? (default: bowtie)
        :type mapper: str -- "bowtie" or "tophat"
        :param index: Which BW index?
        :type index: str -- Pathname to *Bowtie*'s index (just the suffix)
        :param lib_file: Which library file(s)?
        :type lib_file: str
        :param lib_mate_1: Which first mate reads file?
        :type lib_mate_1: str
        :param lib_mate_2: Which second mate reads file?
        :type lib_mate_2: str
        :param qual_file: Which quality file(s)?
        :type qual_file: str
        :param q_mate_1: Which quality file(s) for mate 1?
        :type q_mate_1: str
        :param q_mate_2: Which quality file(s) for mate 2?
        :type q_mate_2: str
        :param annot_file: Which annotation file?
        :type q_mate_2: str
        :param cov_file: Wich coverage file?
        :type cov_file: str
        :param cfg_file: Which configuration file?
        :type cfg_file: str
        :param map_label: A personalized name to the mapping.
        :type map_label: str
        :param barrier: When the pipeline will stop?
        :type barrier: int
        :param quiet: If True don't print anything besides errors.
        :type quiet: bool
        :param id_: Identification number
        :type id_: int

        """
        self.out_dir = out_dir
        self.mapper = mapper
        self.index = index
        self.lib_file = lib_file
        self.lib_mate_1 = lib_mate_1
        self.lib_mate_2 = lib_mate_2
        self.qual_file = qual_file
        self.q_mate_1 = q_mate_1
        self.q_mate_2 = q_mate_2
        self.annot_file = annot_file
        self.cov_file = cov_file
        self.cfg_file =  cfg_file
        self.wedr_barrier = barrier
        self.quiet = quiet
        self.id_ = id_
        self._mapper_out = None or map_label    # Mapper output names.
        self.aln_file = None                    # Alignment file.
        self._mapper_cmd = None                 # Command line of the mapping software.
        self._bb_cmd = None                     # Command line of bowtie-build.
        self._out_pref = None                   # Output prefix.
        self._map_mode = 0                      # Pair or single-end mapping? Use quality files?
        self._bowtie_rrm = 0                    # Bowtie's reads report mode (--al, --un, --max)
        self.log_dir = None

        self.parse_args(out_dir=out_dir,
                        mapper=mapper,
                        index=index,
                        lib_file=lib_file,
                        lib_mate_1=lib_mate_1,
                        lib_mate_2=lib_mate_2,
                        qual_file=qual_file,
                        q_mate_1=q_mate_1,
                        q_mate_2=q_mate_2,
                        annot_file=annot_file,
                        cov_file=cov_file,
                        cfg_file=cfg_file,
                        map_label=map_label,
                        barrier=barrier,
                        quiet=quiet)

    # TODO Arrumar '_parse_args' para lidar com '--just-map', '--just-de' etc.
    def parse_args(self, out_dir="./wedr_out", mapper="bowtie", index=None,
                   lib_file=None, lib_mate_1=None, lib_mate_2=None,
                   qual_file=None, q_mate_1=None, q_mate_2=None, annot_file=None,
                   cov_file=None, cfg_file=None, map_label=None, quiet=False,
                   barrier=0):
        """This method verifies if **Wedring** parameters are set correctly and
        makes other adjustments.

        :param out_dir: Which output directory? (default: ./wedr_out)
        :type out_dir: str -- Pathname to an existing (or not) directory.
        :param mapper: Which mapper? (default: bowtie)
        :type mapper: str -- "bowtie" or "tophat"
        :param index: Which BW index?
        :type index: str -- Pathname to *Bowtie*'s index (just the suffix)
        :param lib_file: Which library file(s)?
        :type lib_file: str
        :param lib_mate_1: Which first mate reads file?
        :type lib_mate_1: str
        :param lib_mate_2: Which second mate reads file?
        :type lib_mate_2: str
        :param qual_file: Which quality file(s)?
        :type qual_file: str
        :param q_mate_1: Which quality file(s) for mate 1?
        :type q_mate_1: str
        :param q_mate_2: Which quality file(s) for mate 2?
        :type q_mate_2: str
        :param annot_file: Which annotation file?
        :type q_mate_2: str
        :param cfg_file: Which configuration file?
        :type cfg_file: str
        :param map_label: A personalized name to the mapping.
        :type map_label: str
        :param barrier: When the pipeline will stop?
        :type barrier: int
        :param quiet: If True don't print anything besides errors.
        :type quiet: bool
        :raises: :class:WedringError

        """
        if self.wedr_barrier != 3:
            mp = mapper.lower()
            if mp in ("bowtie", "tophat"):
                if wedr_check_program(mp):
                    self.mapper = mp
            else:
                raise WedringError(134, "Invalid mapper name: %s" % mapper)
            if quiet:
                self.quiet = True
            if out_dir != "./wedr_out":
                self.out_dir = out_dir
                self.log_dir = path.join(out_dir, "log")
            else:
                self.log_dir = path.join(out_dir, "log")
            if cfg_file is not None:
                if wedr_check_path(cfg_file):
                    self.cfg_file = cfg_file
            if index is not None:
                if wedr_check_path(index + ".*"):
                    self.index = index
            else:
                raise WedringError(135, "You must provide Bowtie's BW index.")
            if lib_file is not None:
                lib_temp = [lf for lf in lib_file.split(',') if lf != '']
                if wedr_check_path(lib_temp):
                    self.lib_file = ','.join(lib_temp)
                    self._map_mode = 2
                if qual_file is not None:
                    qual_temp = [qf for qf in qual_file.split(',') if qf != '']
                    if wedr_check_path(qual_temp):
                        if len(lib_temp) != len(qual_temp):
                            raise WedringError(140, "Unbalanced number of library and quality files.")
                        self.qual_file = qual_file
                        self._map_mode += 5
                else:
                    self._map_mode += 6
            elif lib_mate_1 is not None and lib_mate_2 is not None:
                lib_1_temp = [l1 for l1 in lib_mate_1.split(',') if l1 != '']
                lib_2_temp = [l2 for l2 in lib_mate_2.split(',') if l2 != '']
                if wedr_check_path(lib_1_temp + lib_2_temp):
                    if len(lib_1_temp) != len(lib_2_temp):
                        raise WedringError(140, "Unbalanced number of pair-mate libraries.")
                    if q_mate_1 is not None and q_mate_2 is not None:
                        qual_1_temp = [q1 for q1 in q_mate_1.split(',') if q1 != '']
                        qual_2_temp = [q2 for q2 in q_mate_2.split(',') if q2 != '']
                        if wedr_check_path(qual_1_temp + qual_2_temp):
                            if (len(lib_1_temp) != len(qual_1_temp) or
                                len(lib_1_temp) != len(qual_2_temp)):
                                raise WedringError(140, "Unbalanced number of pair-mate libraries and its qualities.")
                            self.lib_mate_1 = ','.join(lib_1_temp)
                            self.lib_mate_2 = ','.join(lib_2_temp)
                            self.q_mate_1 = ','.join(qual_1_temp)
                            self.q_mate_2 = ','.join(qual_2_temp)
                            self._map_mode = 4
                    elif (q_mate_1 is None) ^ (q_mate_2 is None):
                        raise WedringError(135, "You must set both pair-mate quality files.")
                    self.lib_mate_1 = ','.join(lib_1_temp)
                    self.lib_mate_2 = ','.join(lib_2_temp)
                    self._map_mode = 5
            elif (lib_mate_1 is None) ^ (lib_mate_2 is None):
                raise WedringError(135, "You must set both pair-mate library files.")
            else:
                raise WedringError(135, "You must provide the library file(s).")
            if map_label is not None:
                self._mapper_out = path.join(self.out_dir, map_label)
                self._out_pref = map_label
            elif 7 != self._map_mode != 8:
                if ',' in self.lib_mate_1:
                    map_out_pref = (path.split(path.splitext(self.lib_mate_1[:self.lib_mate_1.index(',')])[0])[1],
                                    path.split(path.splitext(self.lib_mate_2[:self.lib_mate_2.index(',')])[0][1]))
                    if map_out_pref[0][:-2] == map_out_pref[1][:-2]:
                        self._mapper_out = (path.join(self.out_dir,
                                                      map_out_pref[0] +
                                                      "_andothers_vs_" +
                                                      path.split(self.index)[1]))
                        self._out_pref = map_out_pref[0]
                    else:
                        self._mapper_out = (path.join(self.out_dir,
                                                      map_out_pref[0][:-2] + '_' +
                                                      map_out_pref[1][:-2] +
                                                      "_andothers_vs_" +
                                                      path.split(self.index)[1]))
                        self._out_pref = (map_out_pref[0][:-2] + '_' +
                                         map_out_pref[1][:-2] + "_andothers_vs_")
                
                else:
                    map_out_pref = (path.split(path.splitext(self.lib_mate_1)[0])[1],
                                    path.split(path.splitext(self.lib_mate_2)[0])[1])
                    if map_out_pref[0][:-2] == map_out_pref[1][:-2]:
                        self._mapper_out = (path.join(self.out_dir,
                                                       map_out_pref[0][:-2] + "_vs_" +
                                                       path.split(self.index)[1]))
                        self._out_pref = map_out_pref[0][:-2]
                    else:
                        self._mapper_out = (path.join(self.out_dir,
                                                      map_out_pref[0][:-2] + '_' +
                                                      map_out_pref[1][:-2] + "_vs_" +
                                                      path.split(self.index)[1]))
                        self._out_pref = (map_out_pref[0][:-2] + '_' +
                                          map_out_pref[1][:-2])
            else:
                if ',' in self.lib_file:
                    map_out_pref = path.split(path.splitext(self.lib_file[:self.lib_file.index(',')])[0])[1]
                    self._mapper_out = (path.join(self.out_dir, map_out_pref +
                                                  "_andothers_vs_" +
                                                  path.split(self.index)[1]))
                    self._out_pref = map_out_pref + "_andothers_vs_"
                else:
                    map_out_pref = path.split(path.splitext(self.lib_file)[0])[1]
                    self._mapper_out = (path.join(self.out_dir, map_out_pref +
                                                  "_vs_" +
                                                  path.split(self.index)[1]))
                    self._out_pref = map_out_pref
            if annot_file is not None:
                if wedr_check_path(annot_file):
                    self.annot_file = annot_file
            else:
                raise WedringError(135, "You must provide the annotation file.")
        else:
            if annot_file is not None:
                if wedr_check_path(annot_file):
                    self.annot_file = annot_file
            else:
                raise WedringError(135, "You must provide the annotation file.")
            if cov_file is not None:
                if wedr_check_path(cov_file):
                    self.cov_file = cov_file
            else:
                raise WedringError(135, "You must provide the coverage file.")

    def parse_mapper_cmd_line(self):
        """Parser of the mapper command line."""
        if self.mapper == "bowtie":
            self.parse_bowtie_cmd_line()
        else:
            self.parse_tophat_cmd_line()

    def parse_bowtie_cmd_line(self):
        """Parser for the *Bowtie* command line."""
        # Bowtie's read report mode retriever.
        rrm_rtv = lambda s, m: {0: ('', '', ''),
                                1: (" --al %s_mapped" % s, '', ''),
                                3: ('', " --un %s_unmapped" % s, ''),
                                5: ('', '', " --max %s_suppressed" % s),
                                4: (" --al %s_mapped" % s,
                                    " --un %s_unmapped" % s, ''),
                                6: (" --al %s_mapped" % s, '',
                                    " --max %s_suppressed" % s),
                                8: ('', " --un %s_unmapped" % s,
                                    " --max %s_suppressed" % s),
                                9: (" --al %s_mapped" % s,
                                    " --un %s_unmapped" % s,
                                    " --max %s_suppressed" % s)
                               } [m]
        mapper_cmd = self.mapper
        if self.cfg_file is not None:
            mapper_cmd += self.parse_config_file(self.mapper)
        al_, un_, max_ = rrm_rtv(self._mapper_out, self._bowtie_rrm)
        if self._map_mode == 7:
            self._mapper_cmd = ("%s -S -Q %s%s%s%s %s %s %s.sam" %
                                (mapper_cmd, self.qual_file,
                                 al_, un_, max_, self.index,
                                 self.lib_file, self._mapper_out))
        elif self._map_mode == 8:
            self._mapper_cmd = ("%s -S%s%s%s %s %s %s.sam" % 
                                (mapper_cmd, al_, un_, max_, 
                                 self.index, self.lib_file, self._mapper_out))
        elif self._map_mode == 4:
            self._mapper_cmd = ("%s -S --Q1 %s --Q2 %s%s%s%s %s -1 %s -2 %s %s.sam" %
                                (mapper_cmd, self.q_mate_1, self.q_mate_2,
                                 al_, un_, max_, self.index, self.lib_mate_1,
                                 self.lib_mate_2, self._mapper_out))
        elif self._map_mode == 5:
            self._mapper_cmd = ("%s -S%s%s%s %s -1 %s -2 %s %s.sam" %
                                (mapper_cmd, al_, un_, max_, self.index,
                                 self.lib_mate_1, self.lib_mate_2,
                                 self._mapper_out))

    def parse_tophat_cmd_line(self):
        """Parser for the *TopHat* command line."""
        mapper_cmd = self.mapper
        if self.cfg_file is not None:
            mapper_cmd += self.parse_config_file(self.mapper)
        if self._map_mode == 7:
            self._mapper_cmd = ("%s -Q -o %s %s %s %s" %
                                (mapper_cmd, self._mapper_out, self.index,
                                 self.lib_file, self.qual_file))
        elif self._map_mode == 8:
            self._mapper_cmd = ("%s -o %s %s %s" %
                                (mapper_cmd, self._mapper_out, self.index,
                                 self.lib_file))
        elif self._map_mode == 4:
            self._mapper_cmd = ("%s -Q -o %s %s %s %s %s %s" %
                                (mapper_cmd, self._mapper_out, self.index,
                                 self.lib_mate_1, self.lib_mate_2,
                                 self.q_mate_1, self.q_mate_2))
        elif self._map_mode == 5:
            self._mapper_cmd = ("%s -o %s %s %s %s" %
                                (mapper_cmd, self._mapper_out, self.index,
                                 self.lib_mate_1, self.lib_mate_2))

    def parse_config_file(self, section):
        """This method parses the Wedring configuration file and customizes the
        mapping software (*Bowtie* or *TopHat*) parameters.

        :param section: Section name (*bowtie* or *tophat*)
        :type section: str

        """
        if section == "bowtie":
            return self.parse_bowtie_options()
        elif section == "tophat":
            return self.parse_tophat_options()

    def parse_bowtie_options(self):
        """Parser for the *bowtie* section of the configuration file.

        :returns: Additional options of *Bowtie*'s command line.

        """
        mapper_cmd = ""
        params = {"query_in_file": "fastq",
                  "color": "false",
                  "quals": "none",
                  "Q1": "none",
                  "Q2": "none",
                  "skip": "0",
                  "qupto": "0",
                  "trim5": "0",
                  "trim3": "0",
                  "input_qualities": "phred33",
                  "phred33_quals": "true",
                  "phred64_quals": "false",
                  "solexa_quals": "false",
                  "solexa1.3_quals": "false",
                  "integer_quals": "false",
                  "v": "0",
                  "n": "2",
                  "maqerr": "70",
                  "l": "28",
                  "nomaqround": "false",
                  "minins": "0",
                  "maxins": "250",
                  "mate_orientation": "fr",
                  "not_aln_strand": "none",
                  "maxbts": "125",
                  "pairtries": "100",
                  "tryhard": "false",
                  "chunkmbs": "64",
                  "k": "1",
                  "a": "false",
                  "m_": "0",
                  "_M": "0",
                  "best": "false",
                  "strata": "false",
                  "time": "false",
                  "quiet": "false",
                  "al": "false",
                  "un": "false",
                  "max": "false",
                  "snpphred": "30",
                  "snpfrac": "0.001",
                  "col_keepends": "false",
                  "mapq": "255",
                  "offrate": "0",
                  "threads": "1",
                  "mm": "false",
                  "shmen": "false",
                  "seed": "0"}
        cf_parser = RawConfigParser()
        cf_parser.read(self.cfg_file)
        for param, val in cf_parser.items(self.mapper):
            params[param] = val
        for param, val in params.iteritems():
            if param == "query_in_file" and "fastq" != val == 'fasta':
                    mapper_cmd += " -f"
            elif param == "color" and val == "true":
                mapper_cmd += " -C"
                if params["snpphred"] != "30":
                    mapper_cmd += " --snpphred %s" % params["snpphred"]
                if params["snpfrac"] != "0.001":
                    mapper_cmd += " --snpphred %s" % params["snpfrac"]
                if params["col_keepends"] != "false":
                    mapper_cmd += " --col-keepends"
            elif param == "skip" and val != "0":
                mapper_cmd += " -s %s" % val
            elif param == "qupto" and val != "0":
                mapper_cmd += " -u %s" % val
            elif param == "trim5" and val != "0":
                mapper_cmd += " -5 %s" % val
            elif param == "trim3" and val != "0":
                mapper_cmd += " -3 %s" % val
            elif param == "input_qualties" and val != "phred33":
                mapper_cmd += " --phred64-quals"
            elif param == "solexa_quals" and val != "false":
                mapper_cmd += " --solexa-quals"
            elif param == "solexa1.3_quals" and val != "false":
                mapper_cmd += " --solexa1.3-quals"
            elif param == "integer_quals" and val != "false":
                mapper_cmd += " --integer-quals"
            elif param == "v" and val != "0":
                mapper_cmd += " -v %s" % val
            elif param == "n" and val != "2":
                mapper_cmd += " -n %s" % val
            elif param == "maqerr" and val != "70":
                mapper_cmd += " -e %s" % val
            elif param == "l" and val != "28":
                mapper_cmd += " -l %s" % val
            elif param == "nomaqround" and val != "false":
                mapper_cmd += " --nomaqround"
            elif param == "minins" and val != "0":
                mapper_cmd += " -I %s" % val
            elif param == "maxins" and val != "250":
                mapper_cmd += " -X %s" % val
            elif param == "mate_orientation" and val != "fr":
                if val in ("rf", "ff"):
                    mapper_cmd += " --%s" % val
            elif param == "not_aln_strand" and val != "none":
                if val in ("nofw", "norc"):
                    mapper_cmd += " -%s" % val
            elif param == "maxbts" and val != "125" and val != "800":
                mapper_cmd += " --maxbts %s" % val
            elif param == "pairtries" and val != "100" and (7 != self._map_mode != 8):
                mapper_cmd += " --pairtries %s" % val
            elif param == "tryhard" and val != "false":
                mapper_cmd += " --tryhard"
            elif param == "k" and val != "1":
                mapper_cmd += " -k %s" % val
            elif param == "a" and val != "false":
                mapper_cmd += " -a"
            elif param == "m_" and val != "0":
                mapper_cmd += " -m %s" % val
            elif param == "_M" and val != "0":
                mapper_cmd += " -M %s" % val
            elif param == "best" and val != "false":
                mapper_cmd += " --best"
                if params["strata"] != "false":
                    mapper_cmd += " --strata"
                if params["chunkmbs"] != "64":
                    mapper_cmd += " --chunkmbs %s" % params["chunkmbs"]
            elif param == "time" and val != "false":
                mapper_cmd += " -t"
            elif param == "quiet" and val != "false":
                mapper_cmd += " --quiet"
            elif param == "al" and val != "false":
                self._bowtie_rrm += 1
            elif param == "un" and val != "false":
                self._bowtie_rrm += 3
            elif param == "max" and val != "false":
                self._bowtie_rrm += 5
            elif param == "mapq" and val != "255":
                mapper_cmd += " --mapq %s" % val
            elif param == "offrate" and val != "0":
                mapper_cmd += " -o %s" % val
            elif param == "threads" and val != "1":
                mapper_cmd += " -p %s" % val
            elif param == "mm" and val != "false":
                mapper_cmd += " --mm"
            elif param == "shmen" and val != "false":
                mapper_cmd += " --shmen"
            elif param == "seed" and val != "0":
                mapper_cmd += " --seed %s" % val
        return mapper_cmd

    def parse_tophat_options(self):
        """Parser for the tophat section of the configuration file.

        :returns: Additional options of *TopHat*'s command line

        """
        mapper_cmd = ""
        params = {"bowtie1": "false",
                  "mate_inner_dist": "0",
                  "mate_std_dev": "20",
                  "min_anchor_length": "8",
                  "splice_mismatches": "0",
                  "min_intron_length": "70",
                  "max_intron_length": "500000",
                  "max_insertion_length": "3",
                  "max_insertion_length": "3",
                  "solexa_quals": "false",
                  "solexa_1.3_quals": "false",
                  "color": "false",
                  "num_threads": "1",
                  "integer_quals": "false",
                  "max_multihits": "20",
                  "report_secondary_hits": "false",
                  "report_discordant_pair_alignments": "false",
                  "no_coverage_search": "false",
                  "coverage_search": "false",
                  "microexon_search": "false",
                  "library_type": "fr-unstranded",
                  "n": "2",
                  "genome_read_mismatches": "2",
                  "read_mismatches": "2 ",
                  "bowtie_n": "false",
                  "segment_mismatches": "2",
                  "segment_length": "25",
                  "min_coverage_intron": "50",
                  "max_coverage_intron": "20000",
                  "min_segment_intron": "50",
                  "max_segment_intron": "500000",
                  "keep_tmp": "false",
                  "zpacker": "gzip",
                  "fusion_search": "false",
                  "raw_juncs": "none",
                  "fusion_anchor_length": "20",
                  "fusion_min_dist": "10000000",
                  "fusion_read_mismatches": "2",
                  "fusion_multireads": "2",
                  "fusion_multipairs": "2",
                  "fusion_ignore_chromosomes": "none",
                  "no_novel_juncs": "false",
                  "G": "false",
                  "transcriptome_index": "none",
                  "transcriptome_only": "false",
                  "transcriptome_max_hits": "0",
                  "prefilter_multihits": "false",
                  "insertions": "none",
                  "deletions": "none",
                  "no_novel_indels": "false"}
        cf_parser = RawConfigParser()
        cf_parser.read(self.cfg_file)
        for param, val in cf_parser.items(self.mapper):
            params[param] = val
        for param, val in params.iteritems():
            if param == "bowtie1" and val == "true":
                mapper_cmd += " --bowtie1"
            elif param == "output_dir" and val != "./tophat_out":
                mapper_cmd += " -o %s" % val
            elif param == "mate_inner_dist" and val != "0":
                mapper_cmd += " -r %s" % val
            elif param == "mate_std_dev" and val != "20":
                mapper_cmd += " --mate-std-dev %s" % val
            elif param == "min_anchor_length" and val != "8":
                mapper_cmd += " -a %s" % val
            elif param == "splice_mismatches" and val != "0":
                mapper_cmd += " -m %s" % val
            elif param == "splice_mismatches" and val != "0":
                mapper_cmd += " -m %s" % val
            elif param == "min_intron_length" and val != "70":
                mapper_cmd += " -i %s" % val
            elif param == "max_intron_length" and val != "500000":
                mapper_cmd += " -I %s" % val
            elif param == "max_insertion_length" and val != "3":
                mapper_cmd += " --max-insertion-length %s" % val
            elif param == "max_deletion_length" and val != "3":
                mapper_cmd += " --max-deletion-length %s" % val
            elif param == "solexa_quals" and val == "true":
                mapper_cmd += " --solexa-quals"
            elif param == "solexa1.3_quals" and val == "true":
                mapper_cmd += " --solexa1.3-quals"
            elif param == "color" and val == "true":
                mapper_cmd += " -C"
            elif (param == "integer_quals" and params["color"] != "true" and
                  val == "true"):
                    mapper_cmd += " --integer-quals %s" % val
            elif param == "num_threads" and val != "1":
                mapper_cmd += " -p %s" % val
            elif param == "max_multihits" and  val != "20":
                mapper_cmd += " -g %s" % val
            elif param == "report_secondary_hits" and val == "true":
                mapper_cmd += " --report-secondary-hits"
            elif param == "report_discordant_pair_alignments" and val == "true":
                mapper_cmd += " --report_discordant_pair_alignments"
            elif param == "no_coverage_search" and val == "true":
                mapper_cmd += " --no-coverage-search"
            elif param == "coverage_search"and val == "true":
                mapper_cmd += " --coverage-search"
            elif param == "coverage_search" and val == "true":
                mapper_cmd += " --coverage-search"
            elif param == "microexon_search" and val == "true":
                mapper_cmd += " --microexon-search"
            elif (param == "library_type" and val != "fr-unstranded" and
                 val in ("fr-firststrand","fr-secondstrand")):
                mapper_cmd += " --library-type %s" % val
            elif param == "n" and val != "2":
                mapper_cmd += " -n %s" % val
            elif param == "genome_read_mismatches" and val != "2":
                mapper_cmd += " --genome-read-mismatches %s" % val
            elif param == "read_mismatches" and val != "2":
                mapper_cmd += " --read-mismatches %s" % val
            elif param == "bowtie_n" and val == "true":
                mapper_cmd += " --bowtie-n"
            elif param == "segment_mismatches" and val != "2":
                mapper_cmd += " --segment-mismatches %s" % val
            elif param == "segment_length" and val != "25":
                mapper_cmd += " --segment-length %s" % val
            elif param == "min_coverage_intron" and val != "50":
                mapper_cmd += " --min-coverage-intron %s" % val
            elif param == "max_coverage_intron" and val != "20000":
                mapper_cmd += " --max-coverage-intron %s" % val
            elif param == "min_segment_intron" and val != "50":
                mapper_cmd += " --min-segment-intron %s" % val
            elif param == "max_segment_intron" and val != "500000":
                mapper_cmd += " --min-segment-intron %s" % val
            elif param == "keep_tmp" and val == "true":
                mapper_cmd += " --keep-tmp"
            elif param == "zpacker" and val != "gzip":
                mapper_cmd += " -z %s" % val
            elif param == "fusion_search" and val == "true":
                mapper_cmd += " --fusion-search"
                if params["raw_juncs"] != "none":
                    mapper_cmd += " -j %s" % params["raw_juncs"]
                if params["fusion_anchor_length"] != "20":
                    mapper_cmd += (" --fusion-anchor-length %s" %
                                   params["fusion_anchor_length"])
                if params["fusion_min_dist"] != "10000000":
                    mapper_cmd += (" --fusion-min-dist %s" %
                                   params["fusion_min_dist"])
                if params["fusion_read_mismatches"] != "2":
                    mapper_cmd += (" --fusion-read-mismatches %s" %
                                   params["fusion_read_mismatches"])
                if params["fusion_multireads"] != "2":
                    mapper_cmd += (" --fusion-multireads %s" %
                                   params["fusion_multireads"])
                if params["fusion_multipairs"] != "2":
                    mapper_cmd += (" --fusion-multipairs %s" %
                                   params["fusion_multipairs"])
                if params["fusion_ignore_chromosomes"] != "none":
                    mapper_cmd += (" --fusion-ignore-chromosomes %s" %
                                   params["fusion_ignore_chromosomes"])
            elif param == "raw_juncs" and val != "none":
                mapper_cmd += " -j %s" % val
            elif param == "no_novel_juncs" and val == "true":
                mapper_cmd += " --no-novel-juncs"
            elif param == "G" and val == "true":
                if self.annot_file != None:
                    mapper_cmd += " -G %s" % self.annot_file
                else:
                    wedr_report("Ignoring TopHat's option -G/--GTF.")
            elif param == "transcriptome_index" and val != "none":
                mapper_cmd += " --transcriptome-index %s" % val
            elif param == "transcriptome_only" and val == "true":
                mapper_cmd += " -T"
            elif param == "transcriptome_max_hits" and val != "0":
                mapper_cmd += " -x %s" % val
            elif param == "prefilter_multihits" and val == "true":
                mapper_cmd += " -M"
            elif param == "insertions" and val != "none":
                mapper_cmd += " --insertions %s" % val
            elif param == "deletions" and val != "none":
                mapper_cmd += " --deletions %s" % val
            elif param == "no_novel_indels" and val == "true":
                mapper_cmd += " --no-novel-indels"
        return mapper_cmd

    def exec_mapping_stage(self):
        """Executor of the Wedring pipeline."""
        if self.quiet:
            wedr_check_path(self.index + '.*')
            self.exec_mapping()
            self.aln_file = (path.join(self._mapper_out, "accepted_hits.bam") if
                             self.mapper == "tophat" else
                             self._mapper_out + ".sam")
            wedr_check_path(self.aln_file)
            if wedr_check_program("samtools"):
                self.exec_samtools()
                if wedr_check_program("bedtools"):
                    self.exec_bedtools()
        else:
            wedr_check_path(self.index + '.*')
            self.exec_mapping()
            self.aln_file = (path.join(self._mapper_out, "accepted_hits.bam") if
                             self.mapper == "tophat" else
                             self._mapper_out + ".sam")
            wedr_check_path(self.aln_file)
            if wedr_check_program("samtools"):
                wedr_report("[%s] Processing aligments with SAMtools." %
                            self._out_pref)
                self.exec_samtools()
                wedr_report("[%s] SAMtools - Done!." % self._out_pref)
                if wedr_check_program("bedtools"):
                    wedr_report("[%s] Calculating mapping coverage with BEDTools." %
                                self._out_pref)
                    self.exec_bedtools()
                    wedr_report("[%s] BEDTools - Done!" % self._out_pref)
        return self

    def exec_mapping(self):
        """Executor of the mapping part of the pipeline.

        :raises: :class:WedringError

        """
        if not self.quiet:
            wedr_report("[%s] Mapping reads against reference genome." %
                        self._out_pref)
        self.parse_mapper_cmd_line()
        if not self.quiet:
            wedr_report("[%s] Command line:\n   %s" % (self._out_pref,
                                                       self._mapper_cmd))
        errfile = path.join(self.log_dir, self._out_pref + "_mapping.log")
        mp = BioSoft(command=self._mapper_cmd, errfile=errfile)
        mp.run()
        if 0 != mp.return_code != -1:
            raise WedringError(141, "[%s] %s exitted with status %d. See log file '%s' for more details." %
                               (self._out_pref, mp.program_name, mp.return_code,
                                mp.errfile))
        wedr_clean(mp.errfile)
        if not self.quiet:
            wedr_report("[%s] Mapping - Done!." % self._out_pref)

    def exec_samtools(self):
        """Executor of the SAMtools part of the pipeline.

        :raises: :class:WedringError

        """
        if self.mapper == "bowtie":
            if not self.quiet:
                wedr_report("[%s] Converting SAM file to BAM file." %
                            self._out_pref)
            sam_in_pref = path.splitext(self.aln_file)[0]
            bam_out = sam_in_pref + ".bam"
            errfile = path.join(self.log_dir, self._out_pref + "_view.log")
            st = BioSoft(command="samtools view -bS -o %s %s.sam" %
                         (bam_out, self._mapper_out), errfile=errfile)
            if not self.quiet:
                wedr_report("[%s] Command line:\n    %s" % (self._out_pref,
                                                             st.command))
            st.run()
            if 0 != st.return_code != -1:
                raise WedringError(141, "[%s] %s exitted with status %d. See log file '%s' for more details." %
                                   (self._out_pref, st.program_name,
                                    st.return_code, st.errfile))
            wedr_clean(st.errfile)
            wedr_clean(self.aln_file, force=True)
            self.aln_file = bam_out
            wedr_check_path(self.aln_file)
            if not self.quiet:
                wedr_report("[%s] Sorting BAM file." % self._out_pref)
            errfile = path.join(self.log_dir, self._out_pref + "_sort.log")
            st = BioSoft(command = "samtools sort %s %s" %
                         (bam_out, sam_in_pref), errfile =errfile)
            if not self.quiet:
                wedr_report("[%s] Command line:\n    %s" % (self._out_pref,
                                                            st.command))
            st.run()
            if 0 != st.return_code != -1:
                raise WedringError(141, "[%s] %s exitted with status %d. See log file '%s' for more details." %
                                   (self._out_pref, st.program_name,
                                    st.return_code, st.errfile))
            wedr_clean(st.errfile)
            wedr_check_path(self.aln_file)
        if not self.quiet:
            wedr_report("[%s] Indexing BAM file." % self._out_pref)
        errfile = path.join(self.log_dir, self._out_pref + "_index.log")
        st = BioSoft(command="samtools index %s" % self.aln_file,
                      errfile=errfile)
        if not self.quiet:
            wedr_report("[%s] Command line:\n    %s" % (self._out_pref,
                                                         st.command))
        st.run()
        if 0 != st.return_code != -1:
            raise WedringError(141, "[%s] %s exitted with status %d. See log file '%s' for more details." %
                               (self._out_pref, st.program_name, st.return_code,
                                st.errfile))
        wedr_clean(st.errfile)
        wedr_check_path(self.aln_file + ".bai")

    def exec_bedtools(self):
        """Executor of the BEDTools part of the pipeline.

        :raises: :class:WedringError

        """
        cov_out = self.aln_file.replace(".bam", ".cov")
        errfile = path.join(self.log_dir, self._out_pref + "_coverage.log")
        bt = BioSoft(command="bedtools coverage -s -abam %s -b %s" %
                     (self.aln_file, self.annot_file),
                     outfile=cov_out, errfile=errfile)
        if not self.quiet:
            wedr_report("[%s] Command line:\n    %s" % (self._out_pref,
                                                         bt.command))
        bt.run()
        if 0 != bt.return_code != -1:
            raise WedringError(141, "[%s] %s exitted with status %d. See log file '%s' for more details." %
                               (self._out_pref, bt.program_name, bt.return_code,
                                bt.errfile))
        wedr_clean(bt.errfile)
        if wedr_check_path(cov_out):
            self.cov_file = cov_out
