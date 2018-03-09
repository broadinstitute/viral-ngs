'''Automated recording of data provenance and metrics.

This module enables automated recording of metadata and metrics in a uniform format into a centralized place.  The unit of recording
is one command (see cmd.py).  When a command is run, if recording is on, we record: the command; its arguments; info about input and
output files (including hash of file contents and various file metadata); whether the command succeeded or failed; 
info about the code version; info about the runtime environment; the duration of the command.

Eventually, this will enable:
   - answering questions like:
      - how was a given data file made?  (by what command, with what parameters, using what code version)
      - Which workflow versions (parameters, code versions) tend to produce the best results for which kinds of input, 
      according to given metrics?
   - avoiding redundant computation, when a command is re-run with the same inputs and same implementation of the command
'''

from ._metadata.recording import (
    set_run_id,
    cmd_handle_file_arg, cmd_configure_parser, cmd_call_cmd
)
