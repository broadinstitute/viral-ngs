# Unit tests for util.file.py

__author__ = "ilya@broadinstitute.org"

import os, os.path
import util.file

def testTempFiles():
    '''Test creation of tempfiles, as well as dump_file/slurp_file routines'''
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

                util.file.touch_empty(fn)
                assert os.path.getsize(fn)==0

                tmp_fns.append(fn)

        assert os.path.isfile(my_tmp_fn) and not os.path.isfile(my_tmp_fns[0]) and not os.path.isfile(my_tmp_fns[1])

    assert not os.path.isfile(my_tmp_fn)
