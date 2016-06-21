#!/usr/bin/env python
''' This script installs every Tool needed by viral-ngs
'''

from __future__ import print_function
import sys
import timeit
import tools
from tools import *

__author__ = "dpark@broadinstitute.org"

def install_all_tools():
    sumtime = 0.0
    n_tools = 0
    n_success = 0
    for tool_class in tools.all_tool_classes():
        t = tool_class()
        print("installing %s .. " % tool_class.__name__, end="", flush=True)
        runtime = timeit.timeit(t.install)
        sumtime += runtime
        success = t.is_installed()
        print("SUCCESS" if success else "FAILED", end="")
        print(" (%0.1f seconds)" % runtime, flush=True)
        if success:
            n_success += 1
        n_tools += 1
    print("Total %d tools attempted, %d succeeded, %d failed, cumulative install time %0.1f seconds" % (
        n_tools, n_success, n_tools - n_success, sumtime))
    return (n_tools == n_success)
        

if __name__ == '__main__':
    sys.exit(1 if install_all_tools() else 0)
