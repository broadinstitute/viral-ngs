'''This gives a number of useful quick methods for dealing with
tab-text files and gzipped files.
'''

__author__ = "dpark@broadinstitute.org"
__version__ = "PLACEHOLDER"
__date__ = "PLACEHOLDER"

import os, gzip, tempfile, shutil, errno, logging
import util.cmd

log = logging.getLogger(__name__)

def get_project_path() :
    '''Return the absolute path of the top-level project, assumed to be the
       parent of the directory containing this script.'''
    # abspath converts relative to absolute path; expanduser interprets ~
    path = __file__                  # path to this script
    path = os.path.expanduser(path)  # interpret ~
    path = os.path.abspath(path)     # convert to absolute path
    path = os.path.dirname(path)     # containing directory: util
    path = os.path.dirname(path)     # containing directory: main project dir
    return path

def get_build_path() :
    '''Return absolute path of "build" directory'''
    return os.path.join(get_project_path(), 'build')

def get_scripts_path() :
    '''Return absolute path of "scripts" directory'''
    return os.path.join(get_project_path(), 'scripts')

def get_test_path() :
    '''Return absolute path of "test" directory'''
    return os.path.join(get_project_path(), 'test')

def get_test_input_path(testClassInstance=None) :
    '''Return the path to the directory for input files for the specified test class'''
    if testClassInstance is not None :
        return os.path.join(get_test_path(), 'input',
                            type(testClassInstance).__name__)
    else:
        return os.path.join(get_test_path(), 'input')

def mkstempfname(suffix='', prefix='tmp', dir=None, text=False):
    ''' There's no other one-liner way to securely ask for a temp file by filename only.
        This calls mkstemp, which does what we want, except that it returns an open
        file handle, which causes huge problems on NFS if we don't close it.  So close
        it first then return the name part only.
    '''
    fd, fn = tempfile.mkstemp(prefix=prefix, suffix=suffix, dir=dir, text=text)
    os.close(fd)
    return fn

def set_tmpDir(name):
    proposed_prefix = ['tmp']
    if name:
        proposed_prefix.append(name)
    for e in ('LSB_JOBID','LSB_JOBINDEX'):
        if e in os.environ:
            proposed_prefix.append(os.environ[e])
    tempfile.tempdir = tempfile.mkdtemp(prefix='-'.join(proposed_prefix)+'-',
                                        dir=util.cmd.find_tmpDir())

def destroy_tmpDir():
    shutil.rmtree(tempfile.tempdir)

def mkdir_p(dirpath):
    ''' Verify that the directory given exists, and if not, create it.
    '''
    try:
        os.makedirs(dirpath)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(dirpath):
            pass
        else: raise

def open_or_gzopen(fname, *opts):
    return fname.endswith('.gz') and gzip.GzipFile(fname, *opts) or open(fname, *opts)

def cat_files_and_header(inFiles, outFile, add_header=None, strip_in_rows=0):
    ''' This is like cat... except that we optionally drop header rows from
        the input files and/or optionally add a new header row to the output
        file.
    '''
    with open_or_gzopen(outFile, 'wt') as outf:
        if add_header:
            outf.write(add_header+'\n')
        for f in inFiles:
            with open_or_gzopen(f, 'rt') as inf:
                for i in range(strip_in_rows):
                    inf.readline()
                for line in inf:
                    outf.write(line)
    return outFile

def intervals(i, n, l):
    ''' Divide something of length l into n equally sized parts and return the
        start-stop values of the i'th part.  Values are 1-based.  Each part
        will be adjacent and non-overlapping with the next part. i must be a
        number from 1 to n.
    '''
    assert 1 <= i <= n and l>=n
    part_size = l//n
    start = 1 + part_size * (i-1)
    stop = part_size * i
    if i==n:
        stop = l
    return (start,stop)

def read_tabfile_dict(inFile):
    ''' Read a tab text file (possibly gzipped) and return contents as an iterator of dicts. '''
    with open_or_gzopen(inFile, 'rt') as inf:
        header = None
        for line in inf:
            row = line.rstrip('\n').split('\t')
            if line.startswith('#'):
                row[0] = row[0][1:]
                header = row
            elif header==None:
                header = row
            else:
                assert len(header)==len(row)
                yield dict((k,v) for k,v in zip(header, row) if v)

def read_tabfile(inFile):
    ''' Read a tab text file (possibly gzipped) and return contents as an iterator of arrays. '''
    with open_or_gzopen(inFile, 'rt') as inf:
        for line in inf:
            if not line.startswith('#'):
                yield line.rstrip('\n').split('\t')

def readFlatFileHeader(filename, headerPrefix='#', delim='\t'):
    with open_or_gzopen(filename, 'rt') as inf:
        header = inf.readline().rstrip('\n').split(delim)
    if header and header[0].startswith(headerPrefix):
        header[0] = header[0][len(headerPrefix):]
    return header

class FlatFileParser:
    ''' Generic flat file parser that parses tabular text input '''
    def __init__(self, lineIter=None, name=None, outType='dict', readHeader=True,
        headerPrefix='#', delim='\t', requireUniqueHeader=False):
        self.lineIter = lineIter
        self.header = None
        self.name = name
        self.headerPrefix = headerPrefix
        self.readHeader = readHeader
        self.delim = delim
        self.requireUniqueHeader = requireUniqueHeader
        self.line_num=0
        assert outType in ('dict','arrayStrict', 'arrayLoose','both')
        self.outType=outType
        assert readHeader or outType in ('arrayStrict', 'arrayLoose')
    def __enter__(self):
        return self
    def __exit__(self, exc_type, exc_val, exc_tb):
        return 0
    def __iter__(self):
        assert self.lineIter
        for row in self.lineIter:
            out = self.parse(row)
            if out != None:
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
        except Exception as e:
            if self.name:
                log.exception("Exception parsing file at line %d.  Line contents: '%s'  File: %s" % (self.line_num, row, self.name))
            else:
                log.exception("Exception parsing file at line %d.  Line contents: '%s'" % (self.line_num, row))
            raise
    def parseHeader(self, row):
        assert row
        self.header = row
        if self.outType != 'arrayLoose':
            assert len(row) == len(dict([(x,0) for x in row]))
    def parseRow(self, row):
        assert self.outType=='arrayLoose' or self.header and len(self.header)==len(row)
        if self.outType=='arrayLoose' or self.outType=='arrayStrict':
            return row
        out = dict([(self.header[i],row[i]) for i in range(len(self.header))])
        if self.outType=='both':
            for i in range(len(self.header)):
                out[i] = row[i]
        return out

class FlatFileMaker:
    def __init__(self, header, mapIter, comment='#', delim='\t', nullchar='', noHeaderline=False):
        self.comment = comment
        self.delim = delim
        self.header = header
        self.nullchar = nullchar
        self.maps = mapIter
        self.noHeaderline = noHeaderline
    def __iter__(self):
        if not self.noHeaderline:
            yield self.comment + self.delim.join(self.header) + "\n"
        for x in self.maps:
            out = [str(x.get(h,self.nullchar)) for h in self.header]
            yield self.delim.join(out) + "\n"


def fastaMaker(seqs, linewidth=60):
    assert linewidth>0
    for id,seq in seqs:
        yield ">"+id+"\n"
        while len(seq)>linewidth:
            line = seq[:linewidth]
            seq = seq[linewidth:]
            yield line+"\n"
        if seq:
            yield seq+"\n"
