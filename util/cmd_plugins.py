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
    """Handle a file arg."""

cmd_plugin_mgr.add_hookspecs(sys.modules[__name__])

#cmd_plugin_mgr.trace.root.setwriter(print)
#undo = cmd_plugin_mgr.enable_tracing()
