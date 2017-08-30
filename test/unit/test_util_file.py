# Unit tests for util.file.py

__author__ = "ilya@broadinstitute.org"

import os, os.path
import pytest
import util.file

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

        


