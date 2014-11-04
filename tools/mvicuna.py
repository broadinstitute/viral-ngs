"tools.Tool for mvicuna."

import os
import tools, util.file

BroadUnixPath = '/gsap/garage-viral/viral/analysis/xyang/programs'\
                '/M-Vicuna/bin/mvicuna'

class MvicunaTool(tools.Tool) :
    def __init__(self, install_methods = None):
        if install_methods == None:
            #path = os.path.join(util.file.get_scripts_path(), 'mvicuna')
            path = BroadUnixPath
            install_methods = [tools.PrexistingUnixCommand(path)]
        tools.Tool.__init__(self, install_methods = install_methods)
