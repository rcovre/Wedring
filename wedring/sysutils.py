"""The **Wedring** pipeline system utilities
This module provides some functions and one class to interact with the operating
system. It is operating system independent.

.. moduleauthor:: Rafael Covre <covrebio06@gmail.com>

"""

from collections import Sequence
from datetime import datetime
from glob import glob
from os import access, environ, mkdir, path, remove, R_OK, W_OK, X_OK
from sys import exit, stderr
from subprocess import call

from wedring.wedrerror import WedringError

__all__ = ["BioSoft", "wedr_check_path", "wedr_check_program", "wedr_quit",
           "wedr_prepare_directory", "wedr_report", "wedr_clean", "wedr_which"]

class BioSoft(object):
    """Class that implements the environment to execute the softwares
    **Wedring** pipeline uses.
    """

    def __init__(self, command, outfile=None, errfile=None):
        """BioSoft constructor.

        :param command: The command to be executed.
        :type command: str
        :param outfile: The name of the output file.
        :type outfile: str
        :param errfile: The name of the error file.
        :type errfile: str

        """
        self.command = command
        self.outfile = outfile
        self.errfile = errfile
        self.program_name = command[: command.find(" ")]
        self.return_code = -1 # -1 means that no call was performed.

    def run(self):
        """Run the command stored at 'command'. Set the variable 'return_code'
        to the return code of the called command =).

        :raises: :class:WedringError

        """
        try:
            out = None if self.outfile is None else open(self.outfile, "w")
            err = None if self.errfile is None else open(self.errfile, "w")

            self.return_code = call(self.command.split(),
                                    stdout=out,
                                    stderr=err)

            if out is not None:
                out.close()
            if err is not None:
                err.close()

        except EnvironmentError as env_err:
            raise WedringError(env_err.errno, "[%s (%d)] %s%s%s." %
                               (type(env_err).__name__,
                                env_err.errno, 
                                env_err.strerror,
                                ': ' if env_err.filename else '',
                                env_err.filename if env_err.filename else ''))

def wedr_quit(message="", status=0, date=True):
    """Exit the current process displaying a message and with a given status.
    Current date is always displayed.

    :param message: The message.
    :type message: str
    :param status: Terminate the process with this value (default 0).
    :type status: int
    :param date: Display time? True or False? (default True)
    :type date: bool

    """
    out_m = message or "Exiting Wedring pipeline."
    wedr_report(out_m, date=date)
    exit(status)

def wedr_report(message="", date=True):
    """Print to stderr a message. Time may also be displayed (format 
    [Www Mmm dd hh:mm:ss yyyy]).

    :param message: The message.
    :type message: str
    :param date: Display time? True or False? (default True)
    :type date: bool

    """
    dt_m = ""
    if date:
        dt_m = "[%s] " % datetime.now().strftime("%c")
    stderr.write("%s%s\n" % (dt_m, message))

def wedr_check_path(pathname):
    """Verify if a path or a list of paths exist and if they refer to one or
    more readable files or directories. In this case return True, otherwise
    raise a 'not found' error. Wildcards may also be used. For example, "Do*"
    will return True if in the current folder there are files or directories
    named "Documents", "Downloads", "Donuts" etc. This strategy is suitable for
    searching Bowtie's EBWT index files.

    :param path: File/directorie(s) path.
    :type path: str
    :raises: :class:WedringError

    """
    if isinstance(pathname, Sequence):
        if isinstance(pathname, str):
            fp_l = glob(pathname)
            if fp_l == []:
                raise WedringError(137, "File '%s' not found." % pathname)
            elif all(path.exists(fp) and access(fp, R_OK) for fp in fp_l):
                return True
            raise WedringError(137, "File '%s' not found or not readable." %
                               pathname)
        else:
            hits = 0
            for pname in pathname:
                fp_l = glob(pname)
                if fp_l == []:
                    raise WedringError(137, "File '%s' not found." % pname)
                elif all(path.exists(fp) and access(fp, R_OK) for fp in fp_l):
                    hits += 1
                else:
                    raise WedringError(137, "File '%s' not found or not readable." %
                                       pname)
            if len(pathname) == hits:
                return True
            raise WedringError(137, "File '%s' not found or not readable." %
                               ','.join(pathname))
    else:
        raise WedringError(134, "Inappropriate file path '{0}'.".format(pathname))

def wedr_check_program(progname):
    """Verify if a path exists and if it refers to an executable file. In this
    case return True, otherwise raise a 'not found' error.

    :param progname:  Program path.
    :type progname: str
    :raises: :class:WedringError

    """
    for path_ in environ["PATH"].split(path.pathsep):
        pbl_path = path.join(path_, progname) # possible path
        if path.isfile(pbl_path) and access(pbl_path, X_OK):
            return True
    if path.isfile(progname) and access(progname, X_OK):
        return True
    raise WedringError(138, "Program '%s' not found or not executable." %
                       progname)

def wedr_prepare_directory(dirname):
    """This function prepares the output directories used by the **Wedring**
    pipeline.

    :param dirname: The name of the directory to be created.
    :type dirname: str
    :raises: :class:WedringError

    """
    if path.isdir(dirname):
        if not access(dirname, W_OK):
            raise WedringError(138, "Directory '%s' already exists, but it's not writable." %
                               dirname)
    else:
        try:
            mkdir(dirname)
        except OSError as os_err:
            raise WedringError(os_err.errno, "[%s (%d)] %s%s%s." %
                               (type(os_err).__name__,
                                os_err.errno,
                                os_err.strerror,
                                ': ' if os_err.filename else '',
                                os_err.filename if os_err.filename else ''))

def wedr_clean(pathname, force=False):
    """This function performs some cleanups during the **Wedring** pipeline
    execution, for example, remove empty files, if it were necessary.

    :param path: Path to a file.
    :type path: str
    :raises: :class:WedringError

    """
    if wedr_check_path(pathname):
        if not force:
            if path.getsize(pathname) == 0:
                try:
                    remove(pathname)
                except OSError as os_err:
                    raise WedringError(os_err.errno, "[%s (%d)] %s%s%s." %
                                       (type(os_err).__name__,
                                       os_err.errno,
                                       os_err.strerror,
                                       ': ' if os_err.filename else '',
                                       os_err.filename if os_err.filename else ''))
        else:
            try:
                remove(pathname)
            except OSError as os_err:
                raise WedringError(os_err.errno, "[%s (%d)] %s%s%s." %
                                   (type(os_err).__name__,
                                   os_err.errno,
                                   os_err.strerror,
                                   ': ' if os_err.filename else '',
                                   os_err.filename if os_err.filename else ''))

def wedr_which(progname):
    """Given a program name, this function returns the absolute path name of the
    program if it exists in the system path.

    :param progname:  Program path.
    :type progname: str
    :returns: Absolute path to program.
    :raises: :class:WedringError

    """
    for env_path in environ["PATH"].split(path.pathsep):
        pbl_path = path.join(env_path, progname) # possible pathname
        if path.exists(pbl_path) and access(pbl_path, R_OK):
            return pbl_path
    raise WedringError(137, "Program '%s' not found or not readable." %
                       progname)
