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

from ._metadata.argparse_arg_types import InFile, OutFile, InFiles, OutFiles, InFilesPrefix, OutFilesPrefix, InFile_OneOf
from ._metadata.recording import (
    # interface with cmd.py
    add_metadata_tracking, add_metadata_arg,
    # for use by pipe/rules/common.rules to create a common ID for a sequence of steps run as part of same pipeline execution
    set_run_id
)
