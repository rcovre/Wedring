"""The **Wedring** Bowtie's Burrows-Wheeler index builder module.

.. moduleauthor:: Rafael Covre <covrebio06@gmail.com>

"""

from os import path
try:
    from ConfigParser import RawConfigParser
except ImportError:
    from configparser import RawConfigParser

from wedring.sysutils import (BioSoft, wedr_check_path, wedr_check_program,
                              wedr_clean, wedr_prepare_directory, wedr_report)
from wedring.wedrerror import WedringError

__all__ = ["IndexBuilder"]

class IndexBuilder(object):
    """This class provides the environment to run *Bowtie*'s Burrows-Wheeler
    index builder, *bowtie-build*.
    """

    def __init__(self, out_dir="./wedr_index", ref_seq=None, index_label=None,
                 cfg_file=None, quiet=False):
        """:class:IndexBuilder constructor.
        
        :param out_dir: Which output directory? (default: ./wedr_out)
        :type out_dir: str
        :param ref_seq: Which reference sequence(s)?
        :type ref_seq: str
        :param index_label: A personalized name to the index.
        :type index_label: str
        :param cfg_file: Which configuration file?
        :type cfg_file: str
        :param quiet: If True don't print anything besides errors.
        :type quiet: bool

        """
        self.out_dir = out_dir
        self.ref_seq = ref_seq
        self.cfg_file = cfg_file
        self.quiet = quiet
        self.index = None
        self._out_pref = None
        self._bb_cmd = None
        self.log_dir = path.join(out_dir, "log")

        self.parse_args(out_dir=out_dir, ref_seq=ref_seq,
                        index_label=index_label, cfg_file=cfg_file, quiet=quiet)

    def parse_args(self, out_dir="./wedr_index", ref_seq=None, index_label=None,
                   cfg_file=None, quiet=False):
        """This method checks the validity of :class:IndexBuilder parameters.


        :param out_dir: Which output directory? (default: ./wedr_out)
        :type out_dir: str
        :param ref_seq: Which reference sequence(s)?
        :type ref_seq: str
        :param index_label: A personalized name to the index.
        :type index_label: str
        :param cfg_file: Which configuration file?
        :type cfg_file: str
        :param quiet: If True don't print anything besides errors.
        :type quiet: bool
        :raises: :class:WedringError 

        """
        if wedr_check_program("bowtie-build"):
            if quiet:
                self.quiet = True
            if out_dir != "./wedr_index":
                self.out_dir = out_dir
            if cfg_file is not None:
                if wedr_check_path(cfg_file):
                    self.cfg_file = cfg_file
            if ref_seq is not None:
                if wedr_check_path(self.ref_seq.split(',')):
                    self.ref_seq = ref_seq
            else:
                raise WedringError(135, "You must provide the reference sequence file(s).")
            if index_label is not None:
                self.index = path.join(self.out_dir, index_label)
                self._out_pref = index_label
            elif ',' in self.ref_seq:
                index_name = [i for i in self.ref_seq.split(',') if i != ''][0]
                self.index = path.join(self.out_dir,
                                       path.split(path.splitext(index_name)[0])[1] +
                                                    "_andothers")
                self._out_pref = index_name
            else:
                self._out_pref = path.split(path.splitext(self.ref_seq)[0])[1]
                self.index = path.join(self.out_dir, self._out_pref)

    def parse_bb_cmd_line(self):
        """This method set **bowtie-build**'s command line.

        """
        bb_cmd = "bowtie-build"
        if self.cfg_file is not None:
            bb_cmd += self.parse_bowtie_build_options()
        self._bb_cmd = "%s %s %s" % (bb_cmd, self.ref_seq, self.index)

    def parse_bowtie_build_options(self):
        """This method parses the specified configuration file and set
        additional options to the **bowtie-build** command line.

        :returns: str -- The options parsed.

        """
        bb_cmd = ""
        params = {"color": "false",
                  "noauto": "false",
                  "packed": "false",
                  "bmax": "4",
                  "dcv": "1024",
                  "nodc": "false",
                  "noref": "false",
                  "justref": "false",
                  "offrate": "5",
                  "ftabchars": "10",
                  "ntoa": "false",
                  "endianness": "little",
                  "seed": "0",
                  "cutoff": "0",
                  "quiet": "false"}
        cf_parser = RawConfigParser()
        cf_parser.read(self.cfg_file)
        for param, val in cf_parser.items("bowtie-build"):
            params[param] = val
        for param, val in params.iteritems():
            if param == "color" and val != "false":
                bb_cmd += " -C"
            elif param == "noauto" and val != "false":
                bb_cmd += " -a"
                if params["packed"] != "false":
                    bb_cmd += " -p"
                if params["bmax"] != "4":
                    bb_cmd += " --bmax %s" % params["bmax"]
                if params["dcv"] != "1024":
                    bb_cmd += " --dcv %s" % params["dcv"]
            elif param == "nodc" and val != "false":
                bb_cmd += " --nodc"
            elif param == "noref" and val != "false":
                bb_cmd += " -r"
            elif param == "justref" and val != "false":
                bb_cmd += " -3"
            elif param == "offrate" and val != "5":
                bb_cmd += " -o %s" % val
            elif param == "ftabchars" and val != "10":
                bb_cmd += " -t %s" % val
            elif param == "ntoa" and val != "false":
                bb_cmd += " --ntoa"
            elif param == "endianness" and val != "little":
                bb_cmd += " --big"
            elif param == "seed" and val != "0":
                bb_cmd += " --seed %s" % val
            elif param == "cutoff" and val != "0":
                bb_cmd += " --cutoff %s" % val
            elif param == "quiet" and val != "false":
                bb_cmd += " -q"
        return bb_cmd

    def run(self):
        """Execute *bowtie-build*.

        :raises: :class:WedringError

        """
        if not self.quiet:
            wedr_report("[%s] Building Bowtie BW index '%s'." % (self._out_pref,
                                                                  self.index))
        self.parse_bb_cmd_line()
        if not self.quiet:
            wedr_report("[%s] Command line:\n    %s" % (self._out_pref,
                                                         self._bb_cmd))
        wedr_prepare_directory(self.out_dir)
        wedr_prepare_directory(self.log_dir)
        outfile = path.join(self.log_dir, self._out_pref + "_build.log")
        bb = BioSoft(command=self._bb_cmd, outfile=outfile)
        bb.run()
        if 0 != bb.return_code != -1:
            raise WedringError(141, "[%s] %s exitted with status %d. See log file '%s' for more details." %
                               (self._out_pref, bb.program_name,
                                bb.return_code, bb.outfile))
        wedr_clean(outfile)
        if not self.quiet:
            wedr_report("[%s] BW index build - Done!." % self._out_pref)
