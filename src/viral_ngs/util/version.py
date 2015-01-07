''' This gets the git version into python-land
'''

__author__ = "dpark@broadinstitute.org"
__version__ = None

import subprocess, os.path
 
def call_git_describe():
    try:
        cmd = ['git', 'describe', '--tags', '--always', '--dirty']
        out = subprocess.check_output(cmd)
        if type(out) != str:
            out = out.decode('utf-8')
        return out.strip()
    except:
        return None

def release_file():
    path = __file__                  # path to this script
    path = os.path.dirname(path)     # containing directory: util
    path = os.path.dirname(path)     # containing directory: main project dir
    return os.path.join(path, 'VERSION') # the VERSION file
 
def read_release_version():
    try:
        with open(release_file(), 'rt') as inf:
            version = inf.readlines()[0].strip()
    except:
        version = None
    return version
 
def write_release_version(version):
    with open(release_file(), 'wt') as outf:
        outf.write(version+'\n')

def get_version():
    global __version__
    if __version__ == None:
        from_git  = call_git_describe()
        from_file = read_release_version()
        
        if from_git:
            if from_file != from_git:
                write_release_version(from_git)
            __version__ = from_git
        else:
            __version__ = from_file
        
        if __version__ == None:
            raise ValueError("Cannot find the version number!")
    
    return __version__
 
 
if __name__ == "__main__":
    print(get_version())
