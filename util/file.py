'''This gives a number of useful quick methods for dealing with
tab-text files and gzipped files.
'''

__author__ = "dpark@broadinstitute.org"
__version__ = "PLACEHOLDER"
__date__ = "PLACEHOLDER"

import contextlib
import os, os.path
import gzip
import tempfile
import shutil
import errno
import logging
import json
import util.cmd
import util.misc
import sys

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

def check_paths(read=(), write=(), read_and_write=()):
    '''Check that we can read and write the specified files, throw an exception if not.  Useful for checking
    error conditions early in the execution of a command.  Each arg can be a filename or iterable of filenames.
    '''
    read, write, read_and_write = map(util.misc.make_seq, (read, write, read_and_write))
    assert not (set(read) & set(write))
    assert not (set(read) & set(read_and_write))
    assert not (set(write) & set(read_and_write))
    
    for fname in read+read_and_write:
        with open(fname):
            pass

    for fname in write+read_and_write:
        if not os.path.exists(fname):
            with open(fname, 'w'):
                pass
            os.unlink(fname)
        else:
            if not (os.path.isfile(fname) and os.access(fname, os.W_OK)):
                raise PermissionError('Cannot write ' + fname)

def mkstempfname(suffix='', prefix='tmp', directory=None, text=False):
    ''' There's no other one-liner way to securely ask for a temp file by
        filename only.  This calls mkstemp, which does what we want, except
        that it returns an open file handle, which causes huge problems on NFS
        if we don't close it.  So close it first then return the name part only.
    '''
    fd, fn = tempfile.mkstemp(prefix=prefix, suffix=suffix, dir=directory, text=text)
    os.close(fd)
    return fn


@contextlib.contextmanager
def tempfname(*args, **kwargs):
    '''Create a tempfile name on context entry, delete the file (if it exists) on context exit.'''
    fn = mkstempfname(*args, **kwargs)
    try:
        yield fn
    finally:
        if os.path.isfile(fn): os.unlink(fn)


@contextlib.contextmanager
def tempfnames(suffixes, *args, **kwargs):
    '''Create a set of tempfile names on context entry, delete the files (if they exist) on context exit.'''
    fns = [mkstempfname(sfx, *args, **kwargs) for sfx in suffixes]
    try:
        yield fns
    finally:
        if  not keep_tmp():
            for fn in fns:
                if os.path.isfile(fn):
                    os.unlink(fn)

@contextlib.contextmanager
def tmp_dir(*args, **kwargs):
    """Create and return a temporary directory, which is cleaned up on context exit
    unless keep_tmp() is True."""
    try:
        name = tempfile.mkdtemp(*args, **kwargs)
        yield name
    finally:
        if keep_tmp():
            log.debug('keeping tempdir ' + name)
        else:
            shutil.rmtree(name)

@contextlib.contextmanager
def pushd_popd(target_dir):
    '''Temporary change to the specified directory, restoring current directory on context exit.'''
    save_cwd = os.getcwd()
    try:
        os.chdir(target_dir)
        yield target_dir
    finally:
        os.chdir(save_cwd)

def keep_tmp():
    """Whether to preserve temporary directories and files (useful during debugging).
    Return True if the environment variable VIRAL_NGS_TMP_DIRKEEP is set.
    """
    return 'VIRAL_NGS_TMP_DIRKEEP' in os.environ

def set_tmp_dir(name):
    proposed_prefix = ['tmp']
    if name:
        proposed_prefix.append(name)
    for e in ('LSB_JOBID', 'LSB_JOBINDEX', 'JOB_ID'):
        if e in os.environ:
            proposed_prefix.append(os.environ[e])
            break
    tempfile.tempdir = tempfile.mkdtemp(prefix='-'.join(proposed_prefix) + '-', dir=util.cmd.find_tmp_dir())
    os.environ['TMPDIR'] = tempfile.tempdir
    return tempfile.tempdir


def destroy_tmp_dir(tempdir=None):
    if not keep_tmp():
        if tempdir:
            shutil.rmtree(tempdir)
        elif tempfile.tempdir:
            shutil.rmtree(tempfile.tempdir)
    tempfile.tempdir = None


@contextlib.contextmanager
def fifo(num_pipes=1, names=None, name=None):
    pipe_dir = tempfile.mkdtemp()
    pipe_paths = []
    if name is not None:
        names = [name]
    if names:
        num_pipes = len(names)
    for i in range(num_pipes):
        if names is not None:
            fn = names[i]
        else:
            fn = '{}.pipe'.format(i)
        pipe_path = os.path.join(pipe_dir, fn)
        os.mkfifo(pipe_path)
        pipe_paths.append(pipe_path)

    if num_pipes == 1:
        yield pipe_paths[0]
    else:
        yield pipe_paths
    shutil.rmtree(pipe_dir)


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


def open_or_gzopen(fname, *opts, **kwopts):
    # 'U' mode is deprecated in py3 and may be unsupported in future versions, 
    # so use newline=None when 'U' is specified
    for opt in opts:
        if type(opt) == str and 'U' in opt and sys.version_info[0] == 3:
            if 'newline' not in kwopts:
                kwopts['newline'] = None
            break
    return fname.endswith('.gz') and gzip.open(fname, *opts, **kwopts) or open(fname, *opts, **kwopts)


def read_tabfile_dict(inFile):
    ''' Read a tab text file (possibly gzipped) and return contents as an
        iterator of dicts.
    '''
    with open_or_gzopen(inFile, 'rU') as inf:
        header = None
        for line in inf:
            row = [item.strip() for item in line.rstrip('\n').rstrip('\r').split('\t')]
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
    with open_or_gzopen(inFile, 'rU') as inf:
        for line in inf:
            if not line.startswith('#'):
                yield list(item.strip() for item in line.rstrip('\n').rstrip('\r').split('\t'))


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

def make_empty(fname):
    '''Make `fname` a zero-length file with the current timestamp.'''
    with open(fname, 'w'):
        pass

def dump_file(fname, value):
    """store string in file"""
    with open(fname, 'w')  as out:
        out.write(str(value))

def slurp_file(fname, maxSizeMb=50):
    """Read entire file into one string.  If file is gzipped, uncompress it on-the-fly.  If file is larger
    than `maxSizeMb` megabytes, throw an error; this is to encourage proper use of iterators for reading
    large files.  If `maxSizeMb` is None or 0, file size is unlimited."""
    fileSize = os.path.getsize(fname)
    if maxSizeMb  and  fileSize > maxSizeMb*1024*1024:
        raise RuntimeError('Tried to slurp large file {} (size={}); are you sure?  Increase `maxSizeMb` param if yes'.
                           format(fname, fileSize))
    with open_or_gzopen(fname) as f:
        return f.read()

def is_broken_link(filename):
    # isfile() returns True if a file, or a working link
    if os.path.isfile(filename) or os.path.isdir(filename):
        return False
    # otherwise if this is a link
    if os.path.islink(filename):
        # os.path.exists() returns false in the case of broken symlinks
        return not os.path.exists(filename)
    return False


def find_broken_symlinks(rootdir, followlinks=False):
    """
        This function removes broken symlinks within a directory,
        doing the same in each child directory as well (though not following
        functional symlinks, unless they're directories and followlinks=True).
        @param followlinks: only applies to directory links as per os.walk
    """

    broken_links_to_remove = []

    # first check to see if the input is itself a broken link
    if is_broken_link(rootdir):
        broken_links_to_remove.append(rootdir.rstrip("/"))
    else:
        # otherwise traverse the directory hierarchy
        for rootpath, subfolders, files in os.walk(rootdir, followlinks=followlinks):
            for filename in files:
                fpath = os.path.join(rootpath, filename)
                if is_broken_link(fpath):
                    broken_links_to_remove.append(fpath.rstrip("/"))

    return broken_links_to_remove
