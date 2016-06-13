#!/usr/bin/python
''' This gets the git version into python-land
'''

__author__ = "dpark@broadinstitute.org"
__version__ = None

import subprocess
import os
import re
import time, datetime
import os.path
import util.misc


def get_project_path():
    '''Return the absolute path of the top-level project, assumed to be the
       parent of the directory containing this script.'''
    # abspath converts relative to absolute path; expanduser interprets ~
    path = __file__  # path to this script
    path = os.path.expanduser(path)  # interpret ~
    path = os.path.abspath(path)  # convert to absolute path
    path = os.path.dirname(path)  # containing directory: util
    path = os.path.dirname(path)  # containing directory: main project dir
    return path


def call_git_describe():
    cwd = os.getcwd()
    try:
        os.chdir(get_project_path())
        cmd = ['git', 'describe', '--tags', '--always', '--dirty']
        result = util.misc.run_and_print(cmd, silent=True)
        ver = None
        if result.returncode == 0:
            out = result.stdout    
            if not isinstance(out, str):
                out = out.decode('utf-8')
            ver = out.strip()
    except Exception:
        ver = None
    os.chdir(cwd)
    return ver


def release_file():
    return os.path.join(get_project_path(), 'VERSION')


def read_release_version():
    try:
        with open(release_file(), 'rt') as inf:
            version = inf.readlines()[0].strip()
    except Exception:
        version = None
    return version


def write_release_version(version):
    with open(release_file(), 'wt') as outf:
        outf.write(version + '\n')


def approx_version_number():
    """
        In the event that git is unavailable and the VERSION file is not present
        this returns a "version number" in the following precedence:
            - version number from path
                downloads of viral-ngs from GitHub tagged releases
                are likely to be extracted into directories containing
                the version number. If they contain a version number
                in the form d.d.d, we can use it
            - modification time of this file (unix timestamp)
                file modification time for github releases corresponds to
                when the release archives were created, a rough way to ballpark
                the release date. If we can't get the version number from the path
                we can at least use the modification time of this file as a proxy
                for the true version number
            - the current time (unix timestamp)
                the current time is better than not having any version number
    """
    version = ""

    version_re = re.compile(r"(?:(\d+)\.)?(?:(\d+)\.)?(?:(\d+))")
    # path relative to version.py
    viral_ngs_path = os.path.basename(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

    # for tagged releases, it is likely the version number is part of
    # the viral-ngs root directory name
    matches = version_re.search(viral_ngs_path)

    if matches and len([n for n in matches.groups() if n]) == 3:
        version = ".".join( map(str,matches.groups()) )
    else:
        try:
            mtime = os.path.getmtime(__file__)
        except OSError:
            mtime = 0

        if mtime > 0:
            # if we could get the modification time of the current file, use it
            version = str(int(mtime))
        else:
            # just use the current time
            version = str(int(time.time()))

    return version

def get_version():
    global __version__
    if __version__ is None:
        from_git = call_git_describe()
        from_file = read_release_version()

        if from_git:
            if from_file != from_git:
                write_release_version(from_git)
            __version__ = from_git
        else:
            __version__ = from_file

        if __version__ is None:
            __version__ = approx_version_number()
            #raise ValueError("Cannot find the version number!")

    return __version__


if __name__ == "__main__":
    print(get_version())
