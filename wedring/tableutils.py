"""The **Wedring** pipeline table utilities
This module provides some functionalities to generate and write the count
table and the genomic features table to disk. All the functions described here
work on :class:Wedring class instance and/or :class:WedringMast class instance.

.. moduleauthor:: Rafael Covre <covrebio06@gmail.com>

"""

from csv import QUOTE_NONE, reader
from itertools import count, izip, starmap
from operator import eq

from wedring.wedrerror import WedringError

__all__ = ['get_count_table',
           'get_genome_features',
           'prepare_table_header',
           'write_genomic_features_to_file',
           'write_count_table_to_file']

def prepare_table_header(conditions):
    """Prepare the count table header.
    For example, it will turn 'C, C, T, T, T' into 'C1, C2, T1, T2, T3'.

    :param conditions: Comma-separated condition names.
    :type conditions: str
    :returns: A generator which yields strings containing the processed
    header.

    """
    conds_l = [cond for cond in conditions.split(',') if cond != '']
    idx = 0
    conds_l_sz = len(conds_l)
    while idx < conds_l_sz:
        cnt = count(1)
        cur_cond = conds_l[idx]
        while idx < conds_l_sz and conds_l[idx] == cur_cond:
            yield '%s%d' % (conds_l[idx], next(cnt))
            idx += 1

def write_genomic_features_to_file(annotation, filename="genomic_features.txt"):
    """Write to disk a tab-delimited file which contains an identificator and
    the genomic features attributes (gff file ninth column) it corresponds.

    :param annotation: The annotation file.
    :type annotation: str
    :param filename: Name of the resulting file.
    :type filename: str
    :raises: :class:WedringError
    :returns: File name on success.

    """
    try:
        with open(filename, 'w') as ofile:
            ofile.writelines("%s\t%s\n" % (feat[1], feat[5])
                             for feat in get_genome_features(annotation))
    except EnvironmentError as env_err:
        raise WedringError(env_err.errno, "[%s (%d)] %s%s%s." %
                           (type(env_err).__name__,
                            env_err.errno, 
                            env_err.strerror,
                            ': ' if env_err.filename else '',
                            env_err.filename if env_err.filename else ''))
    else:
        return filename

def write_count_table_to_file(annotation, coverage, conditions,
                              filename="count_table.txt"):
    """Write to disk a tab-delimited file corresponding to a count table for
    each genomic feature based on the count information of one or more coverage
    files and the annotation file they were created from.

    :param annotation: The annotation file.
    :type annotation: str
    :param coverage: The coverage file(s).
    :type coverage: 'str' for single file, 'list' or 'tuple' for multiple files.
    :param conditions: Experimental conditions (comma-separated values without
    spaces. Example: 'c,c,t,t,t')
    :type conditions: str
    :param filename: Name of the resulting file.
    :type filename: str
    :raises: :class:WedringError
    :returns: File name on success.

    """
    try:
        with open(filename, 'w') as ofile:
            ofile.writelines('%s\n' % '\t'.join(cnt)
                             for cnt in get_count_table(annotation,
                                                        coverage,
                                                        conditions))
    except EnvironmentError as env_err:
        raise WedringError(env_err.errno, "[%s (%d)] %s%s%s." %
                           (type(env_err).__name__,
                            env_err.errno, 
                            env_err.strerror,
                            ': ' if env_err.filename else '',
                            env_err.filename if env_err.filename else ''))
    else:
        return filename

def get_count_table(annotation, coverage, conditions):
    """Generate the table containing all counts for all genomic features.

    :param annotation: The annotation file.
    :type annotation: str
    :param coverage: The coverage file(s).
    :type coverage: 'str' for single file, 'list' or 'tuple' for multiple
    files.
    :param conditions: Experimental conditions (comma-separated values without
    spaces. Example: 'c,c,t,t,t')
    :type conditions: str
    :returns: A generator which yields a tuple refering to each line of the
    count table.

    """
    feats_lookuptable = hash_genome_features(annotation)
    feats_col = [''] + [feat[1] for feat in get_genome_features(annotation)]
    count_cols = [[header] + ['0'] * (len(feats_col) - 1)
                for header in prepare_table_header(conditions)]
    if isinstance(coverage, list) or isinstance(coverage, tuple):
        cov_list = coverage
    else:
        cov_list = [coverage]
    for column, cov_file in enumerate(cov_list):
        with open(cov_file, 'U') as cov_f:
            for cov_record in reader(cov_f, dialect="excel-tab",
                                     quoting=QUOTE_NONE):
                if '\n' != cov_record[0][0] != '#' and cov_record[9] != '0':
                    for feat in feats_lookuptable[cov_record[8].strip()]:
                        if all(starmap(eq, izip(feat[1:], cov_record[2:5]))):
                            count_cols[column][feat[0]] = cov_record[9]
                            break
    return (i for i in izip(feats_col, *count_cols))

def hash_genome_features(annotation_file):
    """Generate a lookup table for fast searching of the genomic features.

    :param annotation: The annotation file.
    :type annotation: str
    :returns: dict mapped by the genomic feature attribute.

    """
    feats_lookuptable = {}
    _has_key = feats_lookuptable.has_key
    _append = list.append
    for feat in get_genome_features(annotation_file):
        if not _has_key(feat[5]):
            feats_lookuptable[feat[5]] = [feat[:1] + feat[2:5]]
        else:
            _append(feats_lookuptable[feat[5]], feat[:1] + feat[2:5])
    return feats_lookuptable

def get_genome_features(annotation):
    """Get the genomic features from the annotation file.

    :returns: A generator which yields a tuple containing a numeric index
    which the genomic feature corresponds to, the feature, the start and end
    positions and the attributes of the genomic feature.

    .. note::
       For more information about the GFF format see the `GFF
       (General Feature Format) specifications document
       <http://www.sanger.ac.uk/resources/software/gff/spec.htm>`_

    """
    with open(annotation, 'U') as annot_f:
        # Getting maximum number of digits in file lines.
        digits = len(str(sum(1 for _ in annot_f)))
        annot_f.seek(0)
        for idx, feat in enumerate(reader(annot_f, dialect='excel-tab',
                                          quoting=QUOTE_NONE)):
            if '#' != feat[0][0] != '\n':
                feat_no = idx + 1
                feat_id = 'feat%s%d' % ('0' * (digits - len(str(feat_no))),
                                        feat_no)
                feature, start, end, attr = (feat[2], feat[3], feat[4],
                                             feat[8].strip())
                yield feat_no, feat_id, feature, start, end, attr