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

@cmd_hookspec
def cmd_make_arg_handler(mode, compute_fnames):
    """Return an object to pass to the `type` parameter of argparse.ArgumentParser.add_argument(): a callable that takes a string value
    and returns the parsed value."""

cmd_plugin_mgr.add_hookspecs(sys.modules[__name__])
