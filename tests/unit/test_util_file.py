# Unit tests for util.file.py

__author__ = "ilya@broadinstitute.org"

import os, sys, os.path
import builtins
import shutil
import filecmp
import subprocess
import tarfile
import tempfile
import util.file
from test import TestCaseWithTmp
import pytest
from mock import patch

if sys.version_info >= (3, 10, 12):
    '''
    Set the default extraction filter at class instance level for the tarfile
    module, to avoid needing to specify it every time we extract a tarball.
    The extraction filter avoids extracting "unsafe" items from tarballs, and
    is available in Python 3.10.12+. This sets the default to "tar_filter".
    For other filter types and general info on the tarfile module, see:
      https://docs.python.org/3.12/library/tarfile.html#tarfile.tar_filter
      https://docs.python.org/3.12/library/tarfile.html#extraction-filters
    '''
    tarfile.TarFile.extraction_filter = staticmethod(tarfile.tar_filter) if hasattr(tarfile, 'data_filter') else (lambda member, path: member)

def testTempFiles():
    '''Test creation of tempfiles using context managers, as well as dump_file/slurp_file routines'''
    tmp_fns = []
    sfx='tmp-file-test'
    with util.file.tempfname(sfx) as my_tmp_fn:
        tmp_dir_listing = os.listdir( os.path.dirname( my_tmp_fn ) )
        assert os.path.basename( my_tmp_fn ) in tmp_dir_listing

        with util.file.tempfnames((sfx+'-A',sfx+'-B')) as my_tmp_fns:

            assert len(my_tmp_fns)==2
            for fn in my_tmp_fns:
                assert os.path.dirname(fn) == os.path.dirname(my_tmp_fn)
                assert os.path.basename(fn) not in tmp_dir_listing

            for fn in my_tmp_fns:

                assert sfx in fn
                assert os.path.isfile(fn)
                assert os.path.getsize(fn)==0
                assert os.access(fn, os.R_OK | os.W_OK)

                fileValue='my\ntest\ndata\n' + fn + '\n'
                util.file.dump_file( fname=fn, value=fileValue )
                assert os.path.isfile(fn)
                assert os.path.getsize(fn) == len(fileValue)
                assert util.file.slurp_file(fn) == fileValue

                util.file.make_empty(fn)
                assert os.path.getsize(fn)==0

                tmp_fns.append(fn)

        assert os.path.isfile(my_tmp_fn) and not os.path.isfile(my_tmp_fns[0]) and not os.path.isfile(my_tmp_fns[1])

        largeString = 'A' * (2*1024*1024)
        util.file.dump_file(fname=my_tmp_fn, value=largeString)
        with pytest.raises(RuntimeError):
            util.file.slurp_file(my_tmp_fn, maxSizeMb=1)

    assert not os.path.isfile(my_tmp_fn)

def test_check_paths(tmpdir):
    '''Test the util.file.check_paths()'''
    from os.path import join
    from util.file import check_paths
    inDir = util.file.get_test_input_path()
    def test_f(f):
        return join(inDir, f)
    check_paths(read=test_f('empty.bam'))
    check_paths(read=[test_f('empty.bam')])
    with pytest.raises(Exception):
        check_paths(read=test_f('non_existent_file'))
    with pytest.raises(Exception):
        check_paths(write='/non/writable/dir/file.txt')
    writable_dir = str(tmpdir)
    check_paths(write=(join(writable_dir, 'mydata1.txt'),
                       join(writable_dir, 'mydata2.txt')))
    with pytest.raises(Exception):
        check_paths(write=writable_dir)

    util.file.make_empty(join(writable_dir, 'myempty.dat'))
    check_paths(read_and_write=join(writable_dir, 'myempty.dat'))

def test_uncompressed_file_type():
    """Test util.file.uncompressed_file_type()"""
    uft = util.file.uncompressed_file_type
    assert uft('test.fasta.gz') == '.fasta'
    assert uft('test.fasta.bz2') == '.fasta'
    assert uft('test.fasta') == '.fasta'
    assert uft('test.fa.zst') == '.fa'
    assert uft('test.txt.lz4') == '.txt'
    assert uft('test.gz') == ''

    assert uft('/tmp/test.fasta.gz') == '.fasta'
    assert uft('/test/dir/test.fasta.bz2') == '.fasta'
    assert uft('/a/b/c/test.fasta') == '.fasta'
    assert uft('/a/test.gz') == ''

def test_string_to_file_name():
    """Test util.file.string_to_file_name()"""

    unichr = getattr(builtins, 'unichr', chr)

    test_fnames = (
        'simple', 'simple.dat', 'a/b', '/a', '/a/b', '/a/b/', 'a\\b', 
        'a^b&c|d" e::f!g;*?`test`', '(somecmd -f --flag < input > output) && ls',
        'long' * 8000, 'oddchars\\xAA\\o037',
        ''.join(map(chr, range(128))) * 20,
        ''.join(map(unichr, range(1000))),
        )

    with util.file.tmp_dir() as tmp_d:
        for test_fname in test_fnames:
            t_path = os.path.join(tmp_d, util.file.string_to_file_name(test_fname, tmp_d))
            util.file.make_empty(t_path)
            assert os.path.isfile(t_path) and os.path.getsize(t_path) == 0


@pytest.fixture(scope='module', params=['', '.bz2', '.gz', '.lz4', '.zst'])
def compressed_input_file(request):
    return os.path.join(util.file.get_test_input_path(), 'ebola.fasta' + request.param)

@pytest.fixture(scope='module')
def expected_plaintext():
    return os.path.join(util.file.get_test_input_path(), 'ebola.fasta')

def test_decompress_shutil_copyfileobj(request, expected_plaintext, compressed_input_file):
    out = util.file.mkstempfname(os.path.basename(compressed_input_file)+'.fa')
    with util.file.open_or_gzopen(compressed_input_file, 'rt') as inf, open(out, 'wt') as outf:
        shutil.copyfileobj(inf, outf)
    assert filecmp.cmp(out, expected_plaintext, shallow=False)

def test_decompress_line_by_line(request, expected_plaintext, compressed_input_file):
    out = util.file.mkstempfname(os.path.basename(compressed_input_file)+'.fa')
    with util.file.open_or_gzopen(compressed_input_file, 'rt', newline=None) as inf, open(out, 'wt') as outf:
        for line in inf:
            outf.write(line)
    assert filecmp.cmp(out, expected_plaintext, shallow=False)


class TestExtractTarball(TestCaseWithTmp):
    def setUp(self):
        super(TestExtractTarball, self).setUp()
        self.input_dir = os.path.join(util.file.get_test_input_path(), 'TestTarballMerger')
        self.input_mixed_files = list(os.path.join(self.input_dir, "mixed-compressed-input", fn)
                for fn in sorted(os.listdir(os.path.join(self.input_dir, "mixed-compressed-input"))))
        self.expected_outputs = list(os.path.join(self.input_dir, "raw-input", fn)
                for fn in sorted(os.listdir(os.path.join(self.input_dir, "raw-input"))))

    def test_simple_extract(self):
        for tarfile, expected in zip(self.input_mixed_files, self.expected_outputs):
            temp_dir = tempfile.mkdtemp()
            util.file.extract_tarball(tarfile, out_dir=temp_dir)
            out_files = os.listdir(temp_dir)
            self.assertEqual(len(out_files), 1)
            self.assertEqualContents(os.path.join(temp_dir, out_files[0]), expected)

class TestTarballMerger(TestCaseWithTmp):
    def setUp(self):
        super(TestTarballMerger, self).setUp()
        self.input_dir = util.file.get_test_input_path(self)
        self.raw_files = ["file{}".format(x) for x in range(1,5)]
        self.input_tgz_files = list(os.path.join(self.input_dir, "compressed-input", fn)
                for fn in sorted(os.listdir(os.path.join(self.input_dir, "compressed-input"))))
        self.input_mixed_files = list(os.path.join(self.input_dir, "mixed-compressed-input", fn)
                for fn in sorted(os.listdir(os.path.join(self.input_dir, "mixed-compressed-input"))))

    def test_simple_merge(self):
        """
            Simple repack test
        """
        temp_dir = tempfile.gettempdir()
        out_tarball_file = os.path.join(temp_dir,"out.tar.gz")

        util.file.repack_tarballs(out_tarball_file,
                                  self.input_mixed_files
                                  )

        tb = tarfile.open(out_tarball_file)
        tb.extractall(path=temp_dir)

        for i in range(len(self.raw_files)):
            inf = os.path.join(self.input_dir,"raw-input",self.raw_files[i])
            outf = os.path.join(temp_dir,self.raw_files[i])

            self.assertEqualContents(inf, outf)

    def test_merge_with_extract(self):
        """
            Test streaming repack with intermediate extraction to disk
        """
        temp_dir = tempfile.gettempdir()
        out_tarball_file = os.path.join(temp_dir,"out.tar.gz")
        out_extracted_path = os.path.join(temp_dir,"extracted")

        util.file.repack_tarballs(out_tarball_file,
                                  self.input_mixed_files,
                                  extract_to_disk_path=out_extracted_path
                                  )

        tb = tarfile.open(out_tarball_file)
        tb.extractall(path=temp_dir)

        # inspect merged
        for i in range(len(self.raw_files)):
            inf = os.path.join(self.input_dir,"raw-input",self.raw_files[i])
            outf = os.path.join(temp_dir,self.raw_files[i])

            self.assertEqualContents(inf, outf)

        # inspect extracted
        for i in range(len(self.raw_files)):
            inf = os.path.join(self.input_dir,"raw-input",self.raw_files[i])
            outf = os.path.join(out_extracted_path,self.raw_files[i])

            self.assertEqualContents(inf, outf)

    def test_merge_with_extract_repack_from_disk(self):
        """
            Test with repack from disk source after extraction
        """
        temp_dir = tempfile.gettempdir()
        out_tarball_file = os.path.join(temp_dir,"out.tar.gz")
        out_extracted_path = os.path.join(temp_dir,"extracted")

        util.file.repack_tarballs(out_tarball_file,
                                  self.input_mixed_files,
                                  extract_to_disk_path=out_extracted_path,
                                  avoid_disk_roundtrip=False
                                  )

        tb = tarfile.open(out_tarball_file)
        tb.extractall(path=temp_dir)

        # inspect merged
        for i in range(len(self.raw_files)):
            inf = os.path.join(self.input_dir,"raw-input",self.raw_files[i])
            outf = os.path.join(temp_dir,self.raw_files[i])

            self.assertEqualContents(inf, outf)

        # inspect extracted
        for i in range(len(self.raw_files)):
            inf = os.path.join(self.input_dir,"raw-input",self.raw_files[i])
            outf = os.path.join(out_extracted_path,self.raw_files[i])

            self.assertEqualContents(inf, outf)

    def test_piped_in_merge(self):
        """
            Test with streamed input
        """
        temp_dir = tempfile.gettempdir()
        out_tarball_file = os.path.join(temp_dir,"out.tar.gz")

        ps = subprocess.Popen("cat {files}".format(files=' '.join(self.input_tgz_files)).split(), stdout=subprocess.PIPE)
        with patch('sys.stdin', ps.stdout):
            util.file.repack_tarballs(out_tarball_file,
                                      ["-"],
                                      pipe_hint_in="gz")
        ps.wait()

        tb = tarfile.open(out_tarball_file)
        tb.extractall(path=temp_dir)

        for i in range(len(self.raw_files)):
            inf = os.path.join(self.input_dir,"raw-input",self.raw_files[i])
            outf = os.path.join(temp_dir,self.raw_files[i])

            self.assertEqualContents(inf, outf)

    @pytest.fixture(autouse=True)
    def capsys(self, capsys):
        self.capsys = capsys

    def test_piped_out_merge(self):
        """
            Test with streamed output
        """
        temp_dir = tempfile.gettempdir()
        out_tarball_file = os.path.join(temp_dir,"out.tar.gz")

        with open(out_tarball_file, "wb", 0) as outf:
            # temporarily disable pytest's capture of sys.stdout
            with self.capsys.disabled():
                with patch('sys.stdout', outf):
                    util.file.repack_tarballs("-",
                                              self.input_mixed_files,
                                              pipe_hint_out="gz")

        tb = tarfile.open(out_tarball_file)
        tb.extractall(path=temp_dir)

        for i in range(len(self.raw_files)):
            inf = os.path.join(self.input_dir,"raw-input",self.raw_files[i])
            outf = os.path.join(temp_dir,self.raw_files[i])

            self.assertEqualContents(inf, outf)

    def test_merge_piped_in_and_out(self):
        """
            Test with streamed input and output
        """
        temp_dir = tempfile.gettempdir()

        out_tarball_file = os.path.join(temp_dir,"out.tar.gz")

        ps = subprocess.Popen("cat {files}".format(files=' '.join(self.input_tgz_files)).split(), stdout=subprocess.PIPE)
        with patch('sys.stdin', ps.stdout):
            with open(out_tarball_file, "wb", 0) as outf:
                # temporarily disable pytest's capture of sys.stdout
                with self.capsys.disabled():
                    with patch('sys.stdout', outf):
                        util.file.repack_tarballs( "-",
                                                 ["-"],
                                                 pipe_hint_out="gz",
                                                 pipe_hint_in="gz")
        ps.wait()

        tb = tarfile.open(out_tarball_file)
        tb.extractall(path=temp_dir)

        for i in range(len(self.raw_files)):
            inf = os.path.join(self.input_dir,"raw-input",self.raw_files[i])
            outf = os.path.join(temp_dir,self.raw_files[i])

            self.assertEqualContents(inf, outf)
