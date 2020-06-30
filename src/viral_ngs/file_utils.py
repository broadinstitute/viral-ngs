#!/usr/bin/env python3
"""
Utilities for dealing with files.
"""

__author__ = "tomkinsc@broadinstitute.org"
__commands__ = []

import argparse
import csv
import logging

import pandas

import util.cmd
import util.file
import util.misc

log = logging.getLogger(__name__)


# ==============================
# ***  merge_tarballs   ***
# ==============================

def merge_tarballs(out_tarball, in_tarballs, threads=None, extract_to_disk_path=None, pipe_hint_in=None, pipe_hint_out=None):
    ''' Merges separate tarballs into one tarball
        data can be piped in and/or out
    '''
    util.file.repack_tarballs(out_tarball, in_tarballs, threads=threads, extract_to_disk_path=extract_to_disk_path, pipe_hint_in=pipe_hint_in, pipe_hint_out=pipe_hint_out)
    return 0
def parser_merge_tarballs(parser=argparse.ArgumentParser()):
    parser.add_argument(
        'out_tarball', 
        help='''output tarball (*.tar.gz|*.tar.lz4|*.tar.bz2|*.tar.zst|-);
                compression is inferred by the file extension.
        Note: if "-" is used, output will be written to stdout and
         --pipeOutHint must be provided to indicate compression type
         when compression type is not gzip (gzip is used by default).
        ''')
    parser.add_argument(
        'in_tarballs', nargs='+',
        help=('input tarballs (*.tar.gz|*.tar.lz4|*.tar.bz2|*.tar.zst)')
    )
    parser.add_argument('--extractToDiskPath',
                        dest="extract_to_disk_path",
                        help='If specified, the tar contents will also be extracted to a local directory.')
    parser.add_argument('--pipeInHint',
                        dest="pipe_hint_in",
                        default="gz",
                        help='If specified, the compression type used is used for piped input.')
    parser.add_argument('--pipeOutHint',
                        dest="pipe_hint_out",
                        default="gz",
                        help='If specified, the compression type used is used for piped output.')
    util.cmd.common_args(parser, (('threads', None), ('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, merge_tarballs, split_args=True)
    return parser
__commands__.append(('merge_tarballs', parser_merge_tarballs))


def parser_rename_fasta_sequences(parser=argparse.ArgumentParser()):
    parser.add_argument("in_fasta", help="input fasta sequences")
    parser.add_argument("out_fasta", help="output (renamed) fasta sequences")
    parser.add_argument("new_name", help="new sequence base name")
    parser.add_argument(
        "--suffix_always",
        help="append numeric index '-1' to <new_name> if only one sequence exists in <input> (default: %(default)s)",
        default=False,
        action="store_true",
        dest="suffix_always"
    )

    util.cmd.common_args(parser, (('tmp_dir', None), ('loglevel', None), ('version', None)))
    util.cmd.attach_main(parser, main_rename_fasta_sequences)
    return parser
def main_rename_fasta_sequences(args):
    ''' Renames the sequences in a fasta file. Behavior modes:
        1. If input file has exactly one sequence and suffix_always is False,
            then the output file's sequence is named new_name.
        2. In all other cases,
            the output file's sequences are named <new_name>-<i> where <i> is an increasing number from 1..<# of sequences>
    '''
    n_seqs = util.file.fasta_length(args.in_fasta)
    with open(args.in_fasta, 'rt') as inf:
      with open(args.out_fasta, 'wt') as outf:
        if (n_seqs == 1) and not args.suffix_always:
          inf.readline()
          outf.write('>' + args.new_name + '\n')
          for line in inf:
            outf.write(line)
        else:
          i = 1
          for line in inf:
            if line.startswith('>'):
              line = ">{}-{}\n".format(args.new_name, i)
              i += 1
            outf.write(line)

    return 0
__commands__.append(('rename_fasta_sequences', parser_rename_fasta_sequences))

## derived cols

class Adder_Table_Map:
    def __init__(self, tab_file):
        self._mapper = {}
        self._default_col = None
        with open(tab_file, 'rt') as inf:
            reader = csv.DictReader(inf, delimiter='\t')
            self._col_name = reader.fieldnames[0]
            self._orig_cols = reader.fieldnames[1:]
            for row in reader:
                if all(v=='*' for k,v in row.items() if k in self._orig_cols):
                    self._default_col = row.get(self._col_name)
                else:
                    k = self._make_key_str(row)
                    v = row.get(self._col_name, '')
                    log.debug("setting {}={} if {}".format(self._col_name, v, k))
                    self._mapper[k] = v
    def _make_key_str(self, row):
        key_str = ':'.join('='.join((k,row.get(k,''))) for k in self._orig_cols)
        return key_str
    def extra_headers(self):
        return (self._col_name,)
    def modify_row(self, row):
        k = self._make_key_str(row)
        v = self._mapper.get(k)
        if v is None and self._default_col:
           v = row.get(self._default_col, '')
        row[self._col_name] = v
        return row

class Adder_Source_Lab_Subset:
    def __init__(self, restrict_string):
        self._prefix = restrict_string.split(';')[0]
        self._restrict_map = dict(kv.split('=') for kv in restrict_string.split(';')[1].split(':'))
    def extra_headers(self):
        return (self._prefix + '_originating_lab', self._prefix + '_submitting_lab')
    def modify_row(self, row):
        if all((row.get(k) == v) for k,v in self._restrict_map.items()):
            row[self._prefix + '_originating_lab'] = row['originating_lab']
            row[self._prefix + '_submitting_lab']  = row['submitting_lab']
        return row

def parser_tsv_derived_cols(parser=argparse.ArgumentParser()):
    parser.add_argument("in_tsv",  type=str, help="input metadata")
    parser.add_argument("out_tsv", type=str, help="output metadata")
    parser.add_argument("--table_map", type=str, nargs='*', help="Mapping tables. Each mapping table is a tsv with a header. The first column is the output column name for this mapping (it will be created or overwritten). The subsequent columns are matching criteria. The value in the first column is written to the output column. The exception is in the case where all match columns are '*' -- in this case, the value in the first column is the column header name to copy over.")
    parser.add_argument("--lab_highlight_loc", type=str, help="This option copies the 'originating_lab' and 'submitting_lab' columns to new ones including a prefix, but only if they match certain criteria. The value of this string must be of the form prefix;col_header=value:col_header=value. For example, 'MA;country=USA:division=Massachusetts' will copy the originating_lab and submitting_lab columns to MA_originating_lab and MA_submitting_lab, but only for those rows where country=USA and division=Massachusetts.")
    util.cmd.common_args(parser, (('loglevel', None), ('version', None)))
    util.cmd.attach_main(parser, tsv_derived_cols, split_args=True)
    return parser
def tsv_derived_cols(in_tsv, out_tsv, table_map=None, lab_highlight_loc=None):
    ''' Modify metadata table to compute derivative columns on the fly and add or replace new columns
    '''
    adders = []
    if table_map:
        for t in table_map:
            adders.append(Adder_Table_Map(t))
    if lab_highlight_loc:
       adders.append(Adder_Source_Lab_Subset(lab_highlight_loc))

    with open(in_tsv, 'rt') as inf:
        reader = csv.DictReader(inf, delimiter='\t')
        out_headers = reader.fieldnames
        for adder in adders:
            out_headers.extend(adder.extra_headers())

        with open(out_tsv, 'wt') as outf:
            writer = csv.DictWriter(outf, out_headers, delimiter='\t')
            writer.writeheader()
            for row in reader:
                for adder in adders:
                    adder.modify_row(row)
                writer.writerow(row)
__commands__.append(('tsv_derived_cols', parser_tsv_derived_cols))


def parser_tsv_join(parser=argparse.ArgumentParser()):
    parser.add_argument("in_tsvs", nargs="+", type=str, help="input tsvs")
    parser.add_argument("out_tsv", type=str, help="output tsv")
    parser.add_argument("--join_id", type=str, required=True, help="column name to join on")
    util.cmd.common_args(parser, (('loglevel', None), ('version', None)))
    util.cmd.attach_main(parser, tsv_join, split_args=True)
    return parser
def tsv_join(in_tsvs, out_tsv, join_id=None):
    ''' full outer join of tables
    '''

    ''' schaluva implementation:
    df_list = list(pandas.read_csv(f, sep="\t") for f in in_tsvs)
    df_concat = pandas.concat(df_list, axis=0, ignore_index=True, sort=False).fillna('NA')
    pandas.DataFrame.to_csv(df_concat, out_tsv, sep='\t', na_rep='NA', index=False)

    but the above causes types to be parsed and re-outputted slightly differently (e.g. ints become floats)
    and importantly it doesn't actually join on a join key

    attempt to reimplement by hand
    '''

    # prep all readers
    readers = list(csv.DictReader(open(fn, 'rt'), delimiter='\t') for fn in in_tsvs)

    # prep the output header
    header = []
    for reader in readers:
        header.extend(reader.fieldnames)
    header = list(util.misc.unique(header))
    if not join_id or join_id not in header:
        raise Exception()

    # merge everything in-memory
    out_ids = []
    out_row_by_id = {}
    for reader in readers:
        for row in reader:
            row_id = row[join_id]
            row_out = out_row_by_id.get(row_id, {})
            for h in header:
                # prefer non-empty values from earlier files in in_tsvs, populate from subsequent files only if missing
                if not row_out.get(h):
                    row_out[h] = row.get(h, '')
            out_row_by_id[row_id] = row_out
            out_ids.append(row_id)
    out_ids = list(util.misc.unique(out_ids))

    # write output
    with open(out_tsv, 'wt') as outf:
        writer = csv.DictWriter(outf, header, delimiter='\t')
        writer.writeheader()
        writer.writerows(out_row_by_id[row_id] for row_id in out_ids)
__commands__.append(('tsv_join', parser_tsv_join))


# =======================
def full_parser():
    return util.cmd.make_parser(__commands__, __doc__)


if __name__ == '__main__':
    util.cmd.main_argparse(__commands__, __doc__)
