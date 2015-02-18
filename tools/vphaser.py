'''
    V-Phaser 2 variant caller
'''

import logging, subprocess, os
import tools, util.file

log = logging.getLogger(__name__)

class Vphaser2Tool(tools.Tool) :
    def __init__(self, install_methods = None) :
        if install_methods == None :
            path = _get_vphaser_path()
            install_methods = [tools.PrexistingUnixCommand(path)]
        tools.Tool.__init__(self, install_methods = install_methods)

    def execute(self, inBam, outDir, numThreads = None) :
        cmd = [self.install_and_get_path(),
               '-i', inBam,
               '-o', outDir
              ]
        cmdStr = ' '.join(cmd)
        envCopy = os.environ.copy()
        if numThreads != None :
            envCopy['OMP_NUM_THREADS'] = str(numThreads)
            cmdStr = 'OMP_NUM_THREADS=%d ' % numThreads + cmdStr
        log.debug(cmdStr)
        
        # Use check_output instead of check_call so that we get error information
        #    if the executable can't run on travis.
        # Also has the effect of suppressing informational messages from vphaser,
        #    which is probably a good thing.
        try :
            subprocess.check_output(cmd, env = envCopy, stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as ex :
            print(ex.output) # Useful in case of no log handler.
            log.error(ex.output)
            raise

def _get_vphaser_path() :
    uname = os.uname()
    if uname[0] == 'Darwin' :
        osName = 'MacOSX'
    elif uname[0] == 'Linux' and uname[4].endswith('64') :
        osName = 'linux64'
    else :
        log.debug('V-Phaser 2 not implemented for OS {}'.format(
            uname[0] + ' ' + uname[4]))
        return ''
    binariesPath = util.file.get_binaries_path()
    return os.path.join(binariesPath, 'V-Phaser-2.0', osName, 'variant_caller')


"""
Process used to get the files in binaries/V-Phaser-2.0:

wget http://www.broadinstitute.org/software/viral/v_phaser_2/v_phaser_2.zip
unzip
Add "#include <limits.h>" to bam_manip.cpp
Modify src/makefile to create makefile.MacOSX and makefile.linux64
Create linux64 and MacOSX subdirectories

On mac, gcc-4.9 and boost were installed using brew.


# CMake
on linux, "use CMake" (perhaps could instead download from http://www.cmake.org/download/)
on mac, "brew install cmake"

# Bamtools (Note: must use same compiler as V-Phaser 2, otherwise link can fail.)
git clone git://github.com/pezmaster31/bamtools.git
cd bamtools
mkdir build
cd build
on linux "cmake .."
on mac "cmake -DCMAKE_CXX_COMPILER=/usr/local/bin/g++-4.9 -DCMAKE_CC_COMPILER=/usr/local/bin/gcc-4.9 .."
make
cd ../..    # back to V-Phaser-2.0 directory

# boost (only on linux; for mac used brew)
wget http://sourceforge.net/projects/boost/files/latest/download?source=files
tar -xzf boost_1_57_0.tar.gz

# make V-Phaser 2
cd src
make -f makefile.linux64 or makefile.MacOSX

# Cleanup
delete all bamtools stuff
delete all boost stuff
"""
