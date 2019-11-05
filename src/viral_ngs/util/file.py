'''This gives a number of useful quick methods for dealing with
tab-text files and gzipped files, as well as fasta files, plus
general file-handling routines.
'''

__author__ = "dpark@broadinstitute.org"

import codecs
import contextlib
import os
import gzip
import io
import tempfile
import subprocess
import shutil
import errno
import logging
import json
import sys
import io
import csv
import inspect
import tarfile
import itertools
import re
import urllib.request

import util.cmd
import util.misc

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO import FastaIO

import pysam

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
    '''Create a tempfile name on context entry, delete the file (if it exists) on context exit.
    The file is kept for debugging purposes if the environment variable VIRAL_NGS_TMP_DIRKEEP is set.
    '''
    fn = mkstempfname(*args, **kwargs)
    try:
        yield fn
    finally:
        if os.path.isfile(fn) and not keep_tmp():
            os.unlink(fn)


@contextlib.contextmanager
def tempfnames(suffixes, *args, **kwargs):
    '''Create a set of tempfile names on context entry, delete the files (if they exist) on context exit.
    The files are kept for debugging purposes if the environment variable VIRAL_NGS_TMP_DIRKEEP is set.
    '''
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

    _args = inspect.getcallargs(tempfile.mkdtemp, *args, **kwargs)
    length_margin = 6
    for pfx_sfx in ('prefix', 'suffix'):
        if _args[pfx_sfx]:
            _args[pfx_sfx] = string_to_file_name(_args[pfx_sfx], file_system_path=_args['dir'], length_margin=length_margin)
            length_margin += len(_args[pfx_sfx].encode('utf-8'))

    name = None
    try:
        name = tempfile.mkdtemp(**_args)
        yield name
    finally:
        if name is not None:
            if keep_tmp():
                log.debug('keeping tempdir ' + name)
            else:
                shutil.rmtree(name, ignore_errors=True)

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


def extract_tarball(tarfile, out_dir=None, threads=None, compression='auto', pipe_hint=None):
    if not (tarfile == '-' or (os.path.exists(tarfile) and not os.path.isdir(tarfile))):
        raise Exception('file does not exist: %s' % tarfile)
    if out_dir is None:
        out_dir = tempfile.mkdtemp(prefix='extract_tarball-')
    else:
        util.file.mkdir_p(out_dir)
    assert compression in ('gz', 'bz2', 'lz4', 'zip', 'zst', 'none', 'auto')
    if compression is 'auto':
        assert tarfile != '-' or pipe_hint, "cannot autodetect on stdin input unless pipe_hint provided"
        # auto-detect compression type based on file name
        if tarfile=='-':
            lower_fname = pipe_hint
        else:
            lower_fname = os.path.basename(tarfile).lower()
        if lower_fname.endswith('.tar'):
            compression = 'none'
        elif lower_fname.endswith('.zip'):
            compression = 'zip'
        elif lower_fname.endswith('.tgz') or lower_fname.endswith('.tar.gz'):
            compression = 'gz'
        elif lower_fname.endswith('.tar.lz4'):
            compression = 'lz4'
        elif lower_fname.endswith('.tar.bz2'):
            compression = 'bz2'
        elif lower_fname.endswith('.tar.zst'):
            compression = 'zst'
        else:
            raise Exception("unsupported file type: %s" % tarfile)

    if compression == 'zip':
        assert tarfile != '-'
        cmd = ['unzip', '-q', tarfile, '-d', out_dir]
        with open(os.devnull, 'w') as fnull:
            subprocess.check_call(cmd, stderr=fnull)
    else:
        if compression == 'gz':
            decompressor = ['pigz', '-dc', '-p', str(util.misc.sanitize_thread_count(threads))]
        elif compression == 'bz2':
            decompressor = ['lbzip2', '-dc', '-n', str(util.misc.sanitize_thread_count(threads))]
        elif compression == 'lz4':
            decompressor = ['lz4', '-d']
        elif compression == 'zst':
            decompressor = ['zstd', '-d']
        elif compression == 'none':
            decompressor = ['cat']
        untar_cmd = ['tar', '-C', out_dir, '-x']
        if os.getuid() == 0:
            # GNU tar behaves differently when run as root vs normal user
            # we want normal user behavior always
            if 'GNU' in subprocess.check_output(['tar', '--version']).decode('UTF-8'):
                untar_cmd.append('--no-same-owner')
        log.debug("cat {} | {} | {}".format(tarfile, ' '.join(decompressor), ' '.join(untar_cmd)))
        with open(os.devnull, 'w') as fnull:
            if tarfile == '-':
                inf = None
            else:
                inf = open(tarfile, 'rb')
            decompress_proc = subprocess.Popen(decompressor,
                stdin=inf, stdout=subprocess.PIPE)
            untar_proc = subprocess.Popen(untar_cmd,
                stdin=decompress_proc.stdout) #, stderr=fnull)
            if untar_proc.wait():
                raise subprocess.CalledProcessError(untar_proc.returncode, untar_cmd)
            if decompress_proc.wait():
                raise subprocess.CalledProcessError(decompress_proc.returncode, decompressor)
            if inf is not None:
                inf.close()
        log.debug("completed unpacking of {} into {}".format(tarfile, out_dir))

    return out_dir


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


def touch_p(path, times=None):
    '''Touch file, making parent directories if they don't exist.'''
    mkdir_p(os.path.dirname(path))
    touch(path, times=times)


def open_or_gzopen(fname, *opts, **kwargs):
    mode = 'r'
    open_opts = list(opts)
    assert type(mode) == str, "open mode must be of type str"

    # 'U' mode is deprecated in py3 and may be unsupported in future versions,
    # so use newline=None when 'U' is specified
    if len(open_opts) > 0:
        mode = open_opts[0]
        if sys.version_info[0] == 3:
            if 'U' in mode:
                if 'newline' not in kwargs:
                    kwargs['newline'] = None
                open_opts[0] = mode.replace("U","")

    # if this is a gzip file
    if fname.endswith('.gz'):
        # if text read mode is desired (by spec or default)
        if ('b' not in mode) and (len(open_opts)==0 or 'r' in mode):
            # if python 2
            if sys.version_info[0] == 2:
                # gzip.open() under py2 does not support universal newlines
                # so we need to wrap it with something that does
                # By ignoring errors in BufferedReader, errors should be handled by TextIoWrapper
                return io.TextIOWrapper(io.BufferedReader(gzip.open(fname)))

        # if 't' for text mode is not explicitly included,
        # replace "U" with "t" since under gzip "rb" is the
        # default and "U" depends on "rt"
        gz_mode = str(mode).replace("U","" if "t" in mode else "t")
        gz_opts = [gz_mode]+list(opts)[1:]
        return gzip.open(fname, *gz_opts, **kwargs)
    else:
        return open(fname, *open_opts, **kwargs)


def read_tabfile_dict(inFile, header_prefix="#", skip_prefix=None, rowcount_limit=None):
    ''' Read a tab text file (possibly gzipped) and return contents as an
        iterator of dicts.
    '''
    with open_or_gzopen(inFile, 'rU') as inf:
        header = None
        lines_read=0
        for line_no,line in enumerate(inf):
            if line_no==0:
                # remove BOM, if present
                line = line.replace('\ufeff','')
            lines_read+=1
            row = [item.strip() for item in line.rstrip('\r\n').split('\t')]
            # skip empty lines/rows
            if len("".join(row)) == 0 or (skip_prefix and line.startswith(skip_prefix)):
                continue
            if line.startswith(header_prefix):
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
                assert len(header) == len(row), "%s != %s" % (len(header), len(row))
                yield dict((k, v) for k, v in zip(header, row) if v)

            if rowcount_limit and lines_read==rowcount_limit:
                break


def read_tabfile(inFile):
    ''' Read a tab text file (possibly gzipped) and return contents as an
        iterator of arrays.
    '''
    with open_or_gzopen(inFile, 'rU') as inf:
        for line_no,line in enumerate(inf):
            if line_no==0:
                # remove BOM, if present
                line = line.replace('\ufeff','')
            if not line.startswith('#'):
                yield list(item.strip() for item in line.rstrip('\r\n').split('\t'))


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


def bam_is_sorted(bam_file_path):
    # Should perhaps be in samtools.py once it moves to pysam
    samfile = pysam.AlignmentFile(bam_file_path, "rb", check_sq=False)
    if "HD" in samfile.header and "SO" in samfile.header["HD"]:
        return samfile.header["HD"]["SO"] in ("coordinate") # also: "queryname"
    else:
        raise KeyError("Could not locate the SO field in the SAM/BAM file header.")


def concat(inputFilePaths, outputFilePath, append=False):
    '''
        This function creates an output file containing the
        lines present in the input file(s), in the order specified
        by the inputFilePaths list.  If `append` is True,
        appends to the output file instead of overwriting it.
    '''
    with open(outputFilePath, 'a' if append else 'w') as outfile:
        for filePath in util.misc.make_seq(inputFilePaths):
            with open(filePath) as infile:
                for line in infile:
                    outfile.write(line)


def download_file(uriToGet, dest, destFileName=None):
    destDir = os.path.realpath(os.path.expanduser(dest))

    req = urllib.request.urlopen(uriToGet)

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

    for line in urllib.request.urlopen(uriToGet):  # .readlines():
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

def _get_pathconf(file_system_path, param_suffix, default):
    """Return a pathconf parameter value for a filesystem.
    """
    param_str = [s for s in os.pathconf_names if s.endswith(param_suffix)]
    if len(param_str) == 1:
        try:
            return os.pathconf(file_system_path, param_str[0])
        except OSError:
            pass
    return default

def max_file_name_length(file_system_path):
    """Return the maximum valid length of a filename (path component) on the given filesystem."""
    return _get_pathconf(file_system_path, '_NAME_MAX', 80)-1

def max_path_length(file_system_path):
    """Return the maximum valid length of a path on the given filesystem."""
    return _get_pathconf(file_system_path, '_PATH_MAX', 255)-1

def sanitize_id_for_sam_rname(string_in):
    #[0-9A-Za-z!#$%&+./:;?@^_|~-]
    # See character set restrictions in SAM/BAM RNAME spec:
    #   https://samtools.github.io/hts-specs/SAMv1.pdf
    #   [0-9A-Za-z!#$%&+./:;?@^_|~-][0-9A-Za-z!#$%&*+./:;=?@^_|~-]*
    # Here we are being conservative and replacing anything disallowed:
    #   [^0-9A-Za-z!#$%&+./:;?@^_|~-]
    disallowed_char_re = re.compile(r'[^0-9A-Za-z!#$%&+./:;?@^_|~-]')
    string_value = disallowed_char_re.sub("_", string_in)

    # condense runs of underscores
    double_underscore_re = re.compile(r'_{2,}')
    string_value = double_underscore_re.sub("_", string_value)

    # ensure all the character removals did not make the name empty
    string_value = string_value or '_'
    print("sanitizing: %s ====> %s  " % (string_in, string_value))
    return string_value

def write_fasta_with_sanitized_ids(fasta_in, out_filepath):
    with open(out_filepath, "w") as handle:
        fasta_out = FastaIO.FastaWriter(handle, wrap=None)
        fasta_out.write_header()
        for record in SeqIO.parse(fasta_in, "fasta"):
            record.id=sanitize_id_for_sam_rname(record.id)
            fasta_out.write_record(record)
    print("out_filepath",out_filepath)
    print("os.path.dirname(out_filepath)",os.path.dirname(out_filepath))
    print("ls -lah")
    for line in subprocess.check_output(["ls","-lah",os.path.dirname(out_filepath)]).decode("utf-8").split("\n"):
        print(line)
    return out_filepath

@contextlib.contextmanager
def fastas_with_sanitized_ids(input_fasta_paths, use_tmp=False):
    """ Returns a list of file paths for fasta files with
        sanitized IDs 
         ( Suitable for Picard; see: https://github.com/samtools/hts-specs/pull/333 )

        input_fasta_paths is a list of file paths to fasta files

        if use_tmp==False, companion fasta files will be created with ".sanitized_ids.fasta" appended
                           in the same location as the input
        if use_tmp==True, temp files will be written instead
    """
    sanitized_fasta_paths=[]
    if use_tmp:
        with tempfnames(["_{inf_name}".format(inf_name=os.path.basename(inf_path)) for inf_path in [input_fasta_paths]]) as temp_fasta_paths:
            for fasta_in, out_filepath in zip([input_fasta_paths], temp_fasta_paths):
                sanitized_fasta_paths.append(write_fasta_with_sanitized_ids(fasta_in, out_filepath))
            yield sanitized_fasta_paths
    else:
        for fasta_in in [input_fasta_paths]:
            in_fasta_basename = os.path.splitext(os.path.basename(fasta_in))[0]
            out_basedir = os.path.realpath(os.path.dirname(fasta_in))
            new_basename = in_fasta_basename
            if new_basename.lower().endswith('.fa'):
                new_basename = new_basename[:-3] + '.sanitized_ids.fa'
            elif new_basename.lower().endswith('.fasta'):
                new_basename = new_basename[:-6] + '.sanitized_ids.fasta'
            else:
                new_basename = new_basename + '.sanitized_ids.fasta'
            out_filepath = os.path.join(out_basedir,new_basename)
            sanitized_fasta_paths.append(write_fasta_with_sanitized_ids(fasta_in, out_filepath))
        yield sanitized_fasta_paths

class TranspositionError(Exception):
    def __init___(self, *args, **kwargs):
        super(TranspositionError, self).__init__(self, *args, **kwargs)

def transposeChromosomeFiles(inputFilenamesList, sampleRelationFile=None, sampleNameListFile=None):
    ''' Input:  a list of FASTA files representing a genome for each sample.
                Each file contains the same number of sequences (chromosomes, segments,
                etc) in the same order.
                If the parameter sampleRelationFile is specified (as a file path),
                a JSON file will be written mapping sample name to sequence position
                in the output.
        Output: a list of FASTA files representing all samples for each
                chromosome/segment for input to a multiple sequence aligner.
                The number of FASTA files corresponds to the number of chromosomes
                in the genome.  Each file contains the same number of samples
                in the same order.  Each output file is a tempfile.
    '''
    outputFilenames = []

    # open all files
    inputFilesList = [util.file.open_or_gzopen(x, 'r') for x in inputFilenamesList]
    # get BioPython iterators for each of the FASTA files specified in the input
    fastaFiles = [SeqIO.parse(x, 'fasta') for x in inputFilesList]

    # write out json file containing relation of
    # sample name to position in output
    if sampleRelationFile:
        with open(os.path.realpath(sampleRelationFile), "w") as outFile:
            # dict mapping sample->index, zero indexed
            sampleIdxMap = dict((os.path.basename(v).replace(".fasta", ""), k)
                                for k, v in enumerate(inputFilenamesList))
            json.dump(sampleIdxMap, outFile, sort_keys=True, indent=4, separators=(',', ': '))

    if sampleNameListFile:
        with open(os.path.realpath(sampleNameListFile), "w") as outFile:
            sampleNameList = [os.path.basename(v).replace(".fasta", "\n") for v in inputFilenamesList]
            outFile.writelines(sampleNameList)

    # for each interleaved record
    for chrRecordList in itertools.zip_longest(*fastaFiles):
        if any(rec is None for rec in chrRecordList):
            raise TranspositionError("input fasta files must all have the same number of sequences")

        outputFilename = util.file.mkstempfname('.fasta')
        outputFilenames.append(outputFilename)
        with open(outputFilename, "w") as outf:
            # write the corresonding records to a new FASTA file
            SeqIO.write(chrRecordList, outf, 'fasta')

    # close all input files
    for x in inputFilesList:
        x.close()

    return outputFilenames

def string_to_file_name(string_value, file_system_path=None, length_margin=0):
    """Constructs a valid file name from a given string, replacing or deleting invalid characters.
    If `file_system_path` is given, makes sure the file name is valid on that file system.
    If `length_margin` is given, ensure a string that long can be added to filename without breaking length limits.
    """
    replacements_dict = {
        "\\": "-", # win directory separator
        "/": "-", # posix directory separator
        os.sep: "-", # directory separator
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

    # comply with file name length limits
    if file_system_path is not None:
        max_len = max(1, max_file_name_length(file_system_path) - length_margin)
        string_value = string_value[:max_len]
        while len(string_value.encode('utf-8')) > max_len:
            string_value = string_value[:-1]

    # ensure all the character removals did not make the name empty
    string_value = string_value or '_'

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

# used by count_and_sort_barcodes
def count_occurrences_in_tsv(filePath, col=0, noise_chr='.', delimiter='\t', include_noise=False):
    file_occurrence_counts = {}
    with open(filePath) as infile:
        for row in csv.reader(infile, delimiter=delimiter):
            if noise_chr not in row[col] or include_noise:
                file_occurrence_counts[row[col]] = file_occurrence_counts.get(row[col], 0) + 1
    return file_occurrence_counts

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


def uncompressed_file_type(fname):
    """Return the original file extension of either a compressed or an uncompressed file."""
    base, ext = os.path.splitext(fname)
    if ext in ('.gz', '.bz2'):
        base, ext = os.path.splitext(base)
    return ext

def repack_tarballs(out_compressed_tarball,
                    input_compressed_tarballs,
                    extract_to_disk_path=None,
                    extract_numeric_owner=False,
                    avoid_disk_roundtrip=True,
                    ignore_zeros=True,
                    pipe_hint_in=None,
                    pipe_hint_out=None,
                    threads=None):
    threads = util.misc.sanitize_thread_count(threads)

    def choose_compressor(filepath, threads=8):
        return_obj = {}
        filepath = filepath.lower()
        if re.search(r'(\.?tgz|\.?gz)$', filepath):
            compressor = 'pigz {threads}'.format(threads="-p "+str(threads) if threads else "").split()
            return_obj["decompress_cmd"] = compressor + ["-dc"]
            return_obj["compress_cmd"] = compressor + ["-c"]
        elif re.search(r'\.?bz2$', filepath):
            compressor = 'lbzip2 {threads}'.format(threads="-n "+str(threads) if threads else "").split()
            return_obj["decompress_cmd"] = compressor + ["-dc"]
            return_obj["compress_cmd"] = compressor + ["-c"]
        elif re.search(r'\.?lz4$', filepath):
            compressor = ['lz4']
            return_obj["decompress_cmd"] = compressor + ["-dc"]
            return_obj["compress_cmd"] = compressor + ["-c"]
        elif re.search(r'\.?zst$', filepath):
            compressor = ['zstd']
            return_obj["decompress_cmd"] = compressor + ["-d"]
            return_obj["compress_cmd"] = compressor + ["-19"]
        elif re.search(r'\.?tar$', filepath):
            compressor = ['cat']
            return_obj["decompress_cmd"] = compressor
            return_obj["compress_cmd"] = compressor
        else:
            raise IOError("An input file of unknown type was provided: %s" % filepath)
        return return_obj

    class FileDiverter(object):
        """
            This reads bytes from a TarInfo file stream, writes them to a disk file
            and returns the buffered bytes as they are read
        """
        def __init__(self, fileinfo, fileobj, written_mirror_file=None, extract_numeric_owner=False):
            self.written_mirror_file = open(written_mirror_file,"wb")
            self.fileinfo = fileinfo
            self.fileobj = fileobj
            self.extract_numeric_owner = extract_numeric_owner

        def __del__(self):
            self.written_mirror_file.close()

            tar_in.chown(self.fileinfo, self.written_mirror_file.name, self.extract_numeric_owner)
            if not self.fileinfo.issym():
                tar_in.chmod(self.fileinfo, self.written_mirror_file.name)
                tar_in.utime(self.fileinfo, self.written_mirror_file.name)

        def read(self, size):
            assert size is not None

            buf = self.fileobj.read(size)
            self.written_mirror_file.write(buf)
            return buf

    if extract_to_disk_path and not os.path.isdir(extract_to_disk_path):
        mkdir_p(extract_to_disk_path)

    if out_compressed_tarball == "-":
        if not pipe_hint_out:
            raise IOError("cannot autodetect compression for stdoud unless pipeOutHint provided")
        compressor = choose_compressor(pipe_hint_out)["compress_cmd"]
        outfile = None
    else:
        compressor =choose_compressor(out_compressed_tarball)["compress_cmd"]
        outfile = open(out_compressed_tarball, "w")

    out_compress_ps = subprocess.Popen(compressor, stdout=sys.stdout if out_compressed_tarball == "-" else outfile, stdin=subprocess.PIPE)

    tar_out = tarfile.open(fileobj=out_compress_ps.stdin, mode="w|")

    for in_compressed_tarball in input_compressed_tarballs:
        if in_compressed_tarball != "-":
            pigz_ps = subprocess.Popen(choose_compressor(in_compressed_tarball)["decompress_cmd"] + [in_compressed_tarball], stdout=subprocess.PIPE)
        else:
            if not pipe_hint_in:
                raise IOError("cannot autodetect compression for stdin unless pipeInHint provided")
            pigz_ps = subprocess.Popen(choose_compressor(pipe_hint_in)["decompress_cmd"] + [in_compressed_tarball], stdout=subprocess.PIPE, stdin=sys.stdin)
        tar_in = tarfile.open(fileobj=pigz_ps.stdout, mode="r|", ignore_zeros=True)

        fileinfo = tar_in.next()
        while fileinfo is not None:
            if extract_to_disk_path:
                target_path = os.path.normpath(os.path.join(extract_to_disk_path, fileinfo.name).rstrip("/"))
                containing_path = os.path.dirname(target_path)
                mkdir_p(containing_path)

                if avoid_disk_roundtrip and fileinfo.isreg():
                    # avoid disk round trip for regular files (don't re-read from disk to write to new tar)
                    fileobj = tar_in.extractfile(fileinfo)
                    tar_out.addfile(fileinfo, fileobj=FileDiverter(fileinfo, fileobj, written_mirror_file=target_path))
                else:
                    # write  to disk, add to new tarball from disk
                    outfile = tar_in.extract(fileinfo, path=extract_to_disk_path)
                    with pushd_popd(extract_to_disk_path):
                        tar_out.add(fileinfo.name)
            else:
                # if we're not extracting to disk, stream directly between tarballs
                fileobj = tar_in.extractfile(fileinfo)
                tar_out.addfile(fileinfo, fileobj=fileobj)

            fileinfo = tar_in.next()
        pigz_ps.wait()
        tar_in.close()
        if pigz_ps.returncode != 0:
            raise subprocess.CalledProcessError(pigz_ps.returncode, "Call error %s" % pigz_ps.returncode)

    tar_out.close()
    out_compress_ps.stdin.close()
    out_compress_ps.wait()
    if out_compress_ps.returncode != 0:
        raise subprocess.CalledProcessError(out_compress_ps.returncode, "Call error %s" % out_compress_ps.returncode)

    if outfile is not None:
        outfile.close()
