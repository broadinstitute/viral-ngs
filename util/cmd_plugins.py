"""Definitions of plugins used by cmd.py"""

import sys

import pluggy

cmd_hookspec = pluggy.HookspecMarker("viral_ngs_cmd")
cmd_hookimpl = pluggy.HookimplMarker("viral_ngs_cmd")
cmd_plugin_mgr = pluggy.PluginManager("viral_ngs_cmd")

@cmd_hookspec
def cmd_configure_parser(parser):
    """Add any parser opts needed by the plugin"""

@cmd_hookspec
def cmd_call_cmd(cmd_main, args, config):
    """Calls the command."""

@cmd_hookspec(firstresult=True)
def cmd_handle_file_arg(val, mode, compute_fnames):
    """Handle a file arg.

    Args:
        val: the original command-line argument denoting input or output file(s)
        mode: 'r' for input files, 'w' for output files.  All file(s) denoted by an arg
           must have the same mode.
        compute_fnames: function that takes `val` and returns a list of denoted files.
    """

cmd_plugin_mgr.add_hookspecs(sys.modules[__name__])

#cmd_plugin_mgr.trace.root.setwriter(print)
#undo = cmd_plugin_mgr.enable_tracing()
