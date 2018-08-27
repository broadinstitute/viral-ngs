# Unit tests for util.file.py

__author__ = "ilya@broadinstitute.org"

import os, os.path
import builtins
import util.file
import pytest

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
    join = os.path.join
    check_paths = util.file.check_paths
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

    with util.file.fifo() as fifo:
        check_paths(read=fifo)
        check_paths(write=fifo)

def test_hash_file():
    '''Test util.file.hash_file()'''
    assert util.file.hash_file(os.path.join(util.file.get_test_input_path(), 
                                            'G5012.3.mini.bam'), 'sha1')=='582fc7212b4ac1a3fa5dedd331225886056b30f7'

def test_is_file_or_pipe():
    """Test is_file_or_pipe()"""
    ifop = util.file.is_file_or_pipe
    assert not ifop('/')
    assert not ifop('/dev/null')
    assert not ifop('/non/existent/file/asdfasdf')
    assert not ifop('.')

    def check_ifop(targ):
        assert ifop(targ)
        for lnk_fn in (os.symlink, os.link):
            with util.file.tempfname() as lnk:
                os.unlink(lnk)
                lnk_fn(targ, lnk)
                assert ifop(lnk)

    with util.file.tempfname() as tf:
        check_ifop(tf)
        
    with util.file.fifo() as tfifo:
        check_ifop(tfifo)

def test_listdir():
    """Test listdir_full_names() and files_or_pipes_in_dir()"""
    lfn = util.file.listdir_full_names
    fopid = util.file.files_or_pipes_in_dir
    with util.file.tmp_dir() as d:
        assert lfn(d) == fopid(d) == []
        f1 = os.path.join(d, 'f1.dat')
        util.file.make_empty(f1)
        assert lfn(d) == fopid(d) == [f1]
        f2 = os.path.join(d, 'f2.dat')
        util.file.make_empty(f2)
        util.file.touch(f1)
        assert lfn(d) == fopid(d) == [f1, f2]
        with util.file.pushd_popd(d):
            assert lfn('.') == fopid('.') == ['./f1.dat', './f2.dat']
            os.mkdir('d1')
            assert lfn('.') == ['./d1', './f1.dat', './f2.dat']
            assert fopid('.') == ['./f1.dat', './f2.dat']
            os.mkfifo('fifo1')
            assert lfn('.') == ['./d1', './f1.dat', './f2.dat', './fifo1']
            assert fopid('.') == ['./f1.dat', './f2.dat', './fifo1']

def test_json_gz_serialization():
    objs = [ None, 0, 1, '', 'test', [], {}, [1], [1,2], {'0':1}, {'0':{'1':2,'3':[4,5], '6':[7,8]}} ]
    for obj in objs:
        assert util.file.from_json_gz(util.file.to_json_gz(obj)) == obj

    def sets_as_lists(x):
        if isinstance(x, set): return list(x)
        return x

    assert util.file.from_json_gz(util.file.to_json_gz([1, set((2,3))], write_obj=sets_as_lists)) == [1, [2, 3]]

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

