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

    def execute(self, inBam, outDir, numThreads = 1) :
        cmd = [self.install_and_get_path(),
               str(numThreads),
               '-i', inBam,
               '-o', outDir
              ]
        log.debug(' '.join(cmd))
        subprocess.check_call(cmd)


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
    return os.path.join(binariesPath, 'V-Phaser-2.0', osName, 'vphaser')


"""
Process used to build v-phaser on linux64:

# V-Phaser 2 source files
Go to viral-ngs/tools/binaries/V-Phaser-2.0.
Files there were obtained from http://www.broadinstitute.org/software/viral/v_phaser_2/v_phaser_2.zip
Then #include <limits.h> was added to bam_manip.cpp, and makefiles were created for max and linux.

# CMake
use CMake (perhaps could instead download from http://www.cmake.org/download/)

# Bamtools
cd linux64
git clone git://github.com/pezmaster31/bamtools.git
cd bamtools
mkdir build
cd build
cmake ..
make
cd ../../..    # back to V-Phaser-2.0 directory

# boost
wget http://sourceforge.net/projects/boost/files/latest/download?source=files
tar -xzf boost_1_57_0.tar.gz

# make V-Phaser 2
# On linux, this uses the shared bamtools libraries. Future: switch to static libraries
cd src
make -f makefile.linux64

# Cleanup
delete all bamtools stuff except the lib dir and the LICENSE file
delete all boost stuff



Process for building v-phaser on mac:

Do the same as on linux64 except for the following:
- instead of "use Cmake", "brew install cmake"
- cd MacOSX instead of linux64
- make -f makefile.MacOSX
- delete all the bamtools stuff; used static linking so no need to keep the libraries
"""
