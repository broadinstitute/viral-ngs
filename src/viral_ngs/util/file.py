'''This gives a number of useful quick methods for dealing with
tab-text files and gzipped files.
'''

__author__ = "dpark@broadinstitute.org"
__version__ = "PLACEHOLDER"
__date__ = "PLACEHOLDER"

import contextlib
import os
import gzip
import tempfile
import shutil
import errno
import logging
import json
import util.cmd
import util.misc

# imports needed for download_file() and webfile_readlines()
import re
# since py3 split up urllib
try:
    from urllib.request import urlopen # pylint: disable=E0611
except ImportError:
    from urllib2 import urlopen

log = logging.getLogger(__name__)


class StringNotFoundException(Exception):
    """When a substring is not found."""
    pass


def get_project_path():
    '''Return the absolute path of the top-level project, assumed to be the
       parent of the directory containing this script.'''
    # abspath converts relative to absolute path; expanduser interprets ~
    path = __file__  # path to this script
    path = os.path.expanduser(path)  # interpret ~
    path = os.path.abspath(path)  # convert to absolute path
    path = os.path.dirname(path)  # containing directory: util
    path = os.path.dirname(path)  # containing directory: main project dir
    return path


def get_build_path():
    '''Return absolute path of "build" directory'''
    return os.path.join(get_project_path(), 'tools', 'build')


def get_scripts_path():
    '''Return absolute path of "scripts" directory'''
    return os.path.join(get_project_path(), 'tools', 'scripts')


def get_binaries_path():
    '''Return absolute path of "binaries" directory'''
    return os.path.join(get_project_path(), 'tools', 'binaries')


def get_test_path():
    '''Return absolute path of "test" directory'''
    return os.path.join(get_project_path(), 'test')


def get_test_input_path(testClassInstance=None):
    '''Return the path to the directory containing input files for the specified
       test class
    '''
    if testClassInstance is not None:
        return os.path.join(get_test_path(), 'input', type(testClassInstance).__name__)
    else:
        return os.path.join(get_test_path(), 'input')


def get_resources():
    ''' Return the project resources dictionary '''
    jsonfile = os.path.join(get_project_path(), 'resources.json')
    with open(jsonfile, 'rt') as inf:
        resources = json.load(inf)
    return resources


def mkstempfname(suffix='', prefix='tmp', directory=None, text=False):
    ''' There's no other one-liner way to securely ask for a temp file by
        filename only.  This calls mkstemp, which does what we want, except
        that it returns an open file handle, which causes huge problems on NFS
        if we don't close it.  So close it first then return the name part only.
    '''
    fd, fn = tempfile.mkstemp(prefix=prefix, suffix=suffix, dir=directory, text=text)
    os.close(fd)
    return fn


def set_tmp_dir(name):
    proposed_prefix = ['tmp']
    if name:
        proposed_prefix.append(name)
    for e in ('LSB_JOBID', 'LSB_JOBINDEX'):
        if e in os.environ:
            proposed_prefix.append(os.environ[e])
    tempfile.tempdir = tempfile.mkdtemp(prefix='-'.join(proposed_prefix) + '-', dir=util.cmd.find_tmp_dir())
    return tempfile.tempdir


def destroy_tmp_dir(tempdir=None):
    if tempdir:
        shutil.rmtree(tempdir)
    elif tempfile.tempdir:
        shutil.rmtree(tempfile.tempdir)
    tempfile.tempdir = None


def mkdir_p(dirpath):
    ''' Verify that the directory given exists, and if not, create it.
    '''
    try:
        os.makedirs(dirpath)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(dirpath):
            pass
        else:
            raise


def open_or_gzopen(fname, *opts):
    return fname.endswith('.gz') and gzip.open(fname, *opts) or open(fname, *opts)


def read_tabfile_dict(inFile):
    ''' Read a tab text file (possibly gzipped) and return contents as an
        iterator of dicts.
    '''
    with open_or_gzopen(inFile, 'rt') as inf:
        header = None
        for line in inf:
            row = [item.strip() for item in line.rstrip('\n').split('\t')]
            if line.startswith('#'):
                row[0] = row[0][1:]
                header = [item for item in row if len(item)]
            elif header is None:
                header = [item for item in row if len(item)]
            else:
                # if a row is longer than the header
                if len(row) > len(header):
                    # truncate the row to the header length, and only include extra items if they are not spaces 
                    # (takes care of the case where the user may enter an extra space at the end of a row)
                    row = row[:len(header)] + [item for item in row[len(header):] if len(item)]
                assert len(header) == len(row)
                yield dict((k, v) for k, v in zip(header, row) if v)


def read_tabfile(inFile):
    ''' Read a tab text file (possibly gzipped) and return contents as an
        iterator of arrays.
    '''
    with open_or_gzopen(inFile, 'rt') as inf:
        for line in inf:
            if not line.startswith('#'):
                yield list(item.strip() for item in line.rstrip('\n').split('\t'))


def readFlatFileHeader(filename, headerPrefix='#', delim='\t'):
    with open_or_gzopen(filename, 'rt') as inf:
        header = inf.readline().rstrip('\n').split(delim)
    if header and header[0].startswith(headerPrefix):
        header[0] = header[0][len(headerPrefix):]
    return header


class FlatFileParser(object):
    ''' Generic flat file parser that parses tabular text input
    '''

    def __init__(self, lineIter=None, name=None, outType='dict',
                 readHeader=True, headerPrefix='#', delim='\t',
                 requireUniqueHeader=False):
        self.lineIter = lineIter
        self.header = None
        self.name = name
        self.headerPrefix = headerPrefix
        self.readHeader = readHeader
        self.delim = delim
        self.requireUniqueHeader = requireUniqueHeader
        self.line_num = 0
        assert outType in ('dict', 'arrayStrict', 'arrayLoose', 'both')
        self.outType = outType
        assert readHeader or outType in ('arrayStrict', 'arrayLoose')

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        return 0

    def __iter__(self):
        assert self.lineIter
        for row in self.lineIter:
            out = self.parse(row)
            if out is not None:
                yield out

    def parse(self, row):
        self.line_num += 1
        try:
            line = row.rstrip('\n').split(self.delim)
            if self.readHeader:
                if self.headerPrefix and row.startswith(self.headerPrefix):
                    line[0] = line[0][len(self.headerPrefix):]
                    assert not (self.requireUniqueHeader and self.header)
                    self.parseHeader(line)
                    return None
                elif not self.header:
                    self.parseHeader(line)
                    return None
                else:
                    return self.parseRow(line)
            else:
                return self.parseRow(line)
        except Exception:
            template = "Exception parsing file at line {}. Line contents: '{}'"
            message = template.format(self.line_num, row)
            if self.name:
                log.exception("%s  File: %s", message, self.name)
            else:
                log.exception(message)
            raise

    def parseHeader(self, row):
        assert row
        self.header = row
        if self.outType != 'arrayLoose':
            assert len(row) == len(dict([(x, 0) for x in row]))

    def parseRow(self, row):
        assert self.outType == 'arrayLoose' or (self.header and len(self.header) == len(row))

        if self.outType == 'arrayLoose' or self.outType == 'arrayStrict':
            return row
        out = {self.header[i]: row[i] for i in range(len(self.header))}
        if self.outType == 'both':
            for i in range(len(self.header)):
                out[i] = row[i]
        return out


def fastaMaker(seqs, linewidth=60):
    assert linewidth > 0

    for idVal, seq in seqs:
        yield ">{}\n".format(idVal)

        while len(seq) > linewidth:
            line = seq[:linewidth]
            seq = seq[linewidth:]
            yield "{}\n".format(line)

        if seq:
            yield seq + "\n"


def makeFastaFile(seqs, outFasta):
    with open(outFasta, 'wt') as outf:
        for line in fastaMaker(seqs):
            outf.write(line)

    return outFasta


def concat(inputFilePaths, outputFilePath):
    '''
        This function creates an output file containing the
        lines present in the input files, in the order specified
        by the inputFilePaths list.
    '''
    with open(outputFilePath, 'w') as outfile:
        for filePath in inputFilePaths:
            with open(filePath) as infile:
                for line in infile:
                    outfile.write(line)


def download_file(uriToGet, dest, destFileName=None):
    destDir = os.path.realpath(os.path.expanduser(dest))

    req = urlopen(uriToGet)

    if not destFileName:
        m = re.search('filename="(?P<filename>.+)"', req.info()['Content-Disposition'])

        if m:
            destFileName = m.group("filename")
        else:
            destFileName = "file"

    destPath = os.path.join(destDir, destFileName)

    with open(destPath, "wb") as outf:
        while True:
            chunk = req.read(1024)
            if not chunk:
                break
            outf.write(chunk)

    return destPath


def webfile_readlines(uriToGet):

    for line in urlopen(uriToGet):  # .readlines():
        cleanedLine = line.decode("utf-8").strip()
        if len(cleanedLine) > 0:
            yield cleanedLine


def replace_in_file(filename, original, new):
    '''Replace the original string with new in file.

    Raises error if the original is not in the file.
    '''
    with open(filename) as f:
        s = f.read()
    if original not in s:
        raise StringNotFoundException("String '%s' not found." % s)
    s = s.replace(original, new)
    with open(filename, 'w') as f:
        f.write(s)


def cat(output_file, input_files):
    '''Cat list of input filenames to output filename.'''
    with open_or_gzopen(output_file, 'wb') as wfd:
        for f in input_files:
            with open_or_gzopen(f, 'rb') as fd:
                shutil.copyfileobj(fd, wfd, 1024*1024*10)


@contextlib.contextmanager
def temp_catted_files(input_files, suffix=None, prefix=None):
    '''Create a temporary file holding catted contents of input_files.'''
    fn = mkstempfname(suffix=suffix, prefix=prefix)
    try:
        cat(fn, input_files)
        yield fn
    finally:
        os.remove(fn)


def string_to_file_name(string_value):
    replacements_dict = {
        "\\": "-", # win directory separator
        "/": "-", # posix directory separator
        "^": "_", # caret
        "&": "_and_", # background
        "\"": "", # double quotes
        r"'": "", # single quotes
        r":": "_", # colon (problem for ntfs)
        r" ": "_", # spaces
        r"|": "-", # shouldn't confuse a vertical bar for a shell pipe
        r"!": ".", # not a bash operator
        r";": ".", # not a terminator
        r"?": "_", # could be mistaken for a wildcard
        r"*": "_", # could be mistaken for a wildcard
        r"`": "_", # no subshells
        r" -": "_-", # could be mistaken for an argument
        r" --": "_--", # could be mistaken for an argument
        r">": "_", # no redirect chars
        r"<": "_", # no redirect chars
        r"(": "__", # special character
        r")": "__", # special character
        r"\\x": "_", # hex char
        r"\\o": "_", # octal char
        #r"\\u": "", # unicode char
        #"": "", # other illegal strings to replace
    }

    # group of ascii control and non-printable characters
    control_chars = ''.join( map(chr, list(range(0,32)) + list(range(127,160)) ) )
    control_char_re = re.compile('[%s]' % re.escape(control_chars))
    string_value = control_char_re.sub("_", string_value)

    # replacements from the dictionary above
    strs_to_replace_re = re.compile(r'|'.join(re.escape(key) for key in replacements_dict.keys()))
    string_value = strs_to_replace_re.sub(lambda x: replacements_dict.get(x.group(), "_"), string_value)

    # condense runs of underscores
    double_underscore_re = re.compile(r'_{2,}')
    string_value = double_underscore_re.sub("_", string_value)

    # condense runs of dashes
    double_dash_re = re.compile(r'-{2,}')
    string_value = double_dash_re.sub("-", string_value)

    # remove leading or trailing periods (no hidden files (*NIX) or missing file extensions (NTFS))
    string_value = string_value.strip(".")

    return string_value

def grep_count(file_path, to_match, additional_flags=None, fixed_mode=True, starts_with=False):
    '''
        This uses grep for fast counting of strings in a file
    '''
    if not os.path.isfile(file_path) or os.path.getsize(file_path)==0:
        return 0

    env = os.environ.copy()
    env['LC_ALL'] = 'C' #use C locale rather than UTF8 for faster grep

    cmd = ["grep"]
    # '-c' returns the match count
    cmd.append("-c")
    if additional_flags:
        cmd.extend(additional_flags)

    # fixed mode cannot be used with starts_with, since it does not match regular expressions
    # only add the fixed_mode flag if we're not using starts_with
    if not starts_with:
        if fixed_mode:
            cmd.append("-F")
        cmd.append(to_match)
    else:
        cmd.append("^"+to_match)

    cmd.append(file_path)

    number_of_seqs = util.misc.run_and_print(cmd, silent=False, check=True, env=env)
    return int(number_of_seqs.stdout.decode("utf-8").rstrip(os.linesep))

def count_str_in_file(in_file, query_str, starts_with=False):
    if not os.path.isfile(in_file) or os.path.getsize(in_file)==0:
        return 0

    if in_file.endswith('.gz'):
        n = 0
        with gzip.open(in_file, 'rt') as inf:
            if starts_with:
                n = sum(1 for line in inf if line.startswith(query_str))
            else:
                n = sum(1 for line in inf if query_str in line)
        return n
    # use grep count for non-gzip files since it seems to be faster than native on Pythons <3.5
    else:
        return grep_count(in_file, query_str, starts_with=starts_with)

def fasta_length(fasta_path):
    '''
        Count number of records in fasta file
    '''
    return count_str_in_file(fasta_path, '>', starts_with=True)

def count_fastq_reads(inFastq):
    '''
        Count number of reads in fastq file
    '''
    n = line_count(inFastq)
    if n % 4 != 0:
        raise Exception("cannot count reads in a fastq with wrapped lines")
    return n // 4
    # unfortunately, both @ and + are potential quality characters and cannot be used in a
    # fastq counting approach....
    #return count_str_in_file(inFastq, '@', starts_with=True)

def line_count(infname):
    n = 0
    with open_or_gzopen(infname, 'rt') as inf:
        for line in inf:
            n += 1
    return n

def touch(fname, times=None):
    with open(fname, 'a'):
        os.utime(fname, times)
