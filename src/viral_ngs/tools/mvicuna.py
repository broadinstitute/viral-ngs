"tools.Tool for mvicuna."

import os
import tools, util.file

# BroadUnixPath = '/gsap/garage-viral/viral/analysis/xyang/programs'\
#                 '/M-Vicuna/bin/mvicuna'

class MvicunaTool(tools.Tool) :
    def __init__(self, install_methods = None):
        if install_methods == None :
            path = _get_mvicuna_path()
            install_methods = [tools.PrexistingUnixCommand(path)]
        tools.Tool.__init__(self, install_methods = install_methods)


def _get_mvicuna_path() :
    uname = os.uname()
    if uname[0] == 'Darwin' :
        osName = 'MacOSX'
    elif uname[0] == 'Linux' and uname[4].endswith('64') :
        osName = 'linux64'
    else :
        log.debug('mvicuna not implemented for OS {}'.format(
            uname[0] + ' ' + uname[4]))
        return ''
    binariesPath = util.file.get_binaries_path()
    return os.path.join(binariesPath, 'mvicuna', osName, 'mvicuna')

"""
Instructions for building mvicuna on Mac OS X Mavericks:
- Install brew
- brew install homebrew/versions/gcc49
- Copy src files and empty bin directory from:
    /gsap/garage-viral/viral/analysis/xyang/programs/M-Vicuna
- Edit makefile in src directory:
    - Set COMPILER to /usr/local/bin/gcc-4.9
    - Add  -lstdc++ to end of compile command
- make (in src directory)
"""