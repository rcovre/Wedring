import csv
import os

__all__ = ['GffFormatError', 'GffValidator', 'write_validated_gff']

class GffFormatError(Exception):
    """Error class to be raised when a gff format error issue happens."""

    def __init__(self, strerror, filename=None):
        self.strerror = strerror
        self.filename = filename
        if self.filename:
            super(GffFormatError, self).__init__(self.strerror, self.filename)
        else:
            super(GffFormatError, self).__init__(self.strerror)

    def __str__(self):
        if self.filename:
            return '%s: \'%s\'' % (self.strerror, self.filename)
        else:
            return self.strerror

class GffValidator(object):
    """Class used to verify the validity of a gff annotation file."""

    def __init__(self, gff_file):
        self.gff_file = gff_file

    def verify_field_number(self):
        """Verify if all gff records have the same length and if this length
        is 9 (standard number of fields according to the gff format
        specifications).

        :returns: bool
        :raises: :class:GffFormatError

        """

        with open(self.gff_file, 'U') as gff_f:
            for line, record in enumerate(gff_f):
                if '#' != record[0] != '\n':
                    if len(record.split('\t')) != 9:
                        raise GffFormatError('invalid number of fields at line'
                                             ' %d' % (line + 1), self.gff_file)
        return True

    def verify_newline_char(self):
        """verify if the newline character present in the gff file is suitable
        to the system.

        :returns: bool
        :raises: :class:GffFormatError

        """
        with open(self.gff_file, 'U') as gff_f:
            for _ in gff_f:
                if isinstance(gff_f.newlines, tuple):
                    raise GffFormatError('gff file newline character differs '
                                         'from the system.')
                else:
                    if gff_f.newlines != os.linesep:
                        raise GffFormatError('gff file newline character '
                                             'differs from the system.')
        return True

    def check_validity(self):
        """Apply all the class verification methods to the gff file.
        Current validity check:
        * Number of fields per record
        * Newline character

        :returns: bool
        :raises: :class:GffFormatError

        """
        return all((self.verify_field_number(), self.verify_newline_char()))

def write_validated_gff(gff_file, gff_out=None):
    """Writes to disk a validated gff file. If the number of fields in the
    records are invalid an exception is raised. Otherwise, if only the newline
    character differs from the system a new file is written with a suitable
    newline character.

    :param gff_file: Annotation file in gff format.
    :type gff_file: str
    :param gff_out: Name of the file where the validated gff will be written.
    :type gff_out: str
    :returns: True if file doesn't need to be validated, else returns the name
    of the new validated gff file after write it to disk.
    :raises: :class:GffFormatError

    """
    validator = GffValidator(gff_file)
    if validator.verify_field_number():
        try:
            validator.verify_newline_char()
        except GffFormatError:
            with open(gff_file, 'U') as gin:
                if not gff_out:
                    fname, ext = os.path.splitext(gff_file)
                    gff_out = '%s_validated%s' % (fname, ext)
                with open(gff_out, 'w') as gout:
                    gff_reader = csv.reader(gin, dialect='excel-tab',
                                            quoting=csv.QUOTE_NONE)
                    gff_records = ('\t'.join(rec.strip() for rec in record)
                                   for record in gff_reader)
                    gout.write('\n'.join(gff_records) + '\n')
        else:
            return True
    return gff_out