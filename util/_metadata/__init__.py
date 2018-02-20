"""Implementation of metadata recording and analysis.  Public API is in util/metadata.py.

Automated recording of data provenance and metrics.

This module enables the following:
   - answering questions like:
      - how was a given data file made?  (by what command, with what parameters, using what code version)
      - Which workflow versions (parameters, code versions) tend to produce the best results for which kinds of input, 
      according to given metrics?
   - avoiding redundant computation, when a command is re-run with the same inputs

../metadata_utils.py provides external command-line interface for querying provenance data.

The unit of recording in this module is one command (see cmd.py), as opposed to one Python function or one whole workflow.
To use this module, import InFile and OutFile from it; then, when defining argparse arguments for a command, in add_argument()
use "type=InFile" for input files or "type=OutFile" for output files. Until provenance tracking is enabled (see below how),
InFile and OutFile are defined to be simply str.  To enable provenance tracking, set the environment variable
VIRAL_NGS_METADATA_PATH to a writable directory.  Then, whenever a command is run, a file is written to that directory
recording the command, all its parameters, the input and output files, and details of the run environment.
Input and output files are identified by a hash of their contents, so that copied/moved/renamed files can be properly identified.
For each file, we also record metadata such as name, size, and modification date.


: for each data file, automatically keeping track of how it was created -- from which inputs, by which 
code version and with what parameters.


Note that only files actually used are listed in the metadata; optional files not specified in a given command invocation are not.


See data_utils.py for utils that create useful reports from provenance data.

Environment variables used:

   VIRAL_NGS_METADATA_PATH: location to which metadata should be recorded.  If this environment variable is not set, metadata recording
       is disabled, and this module has no effect.

Implementation notes:

Metadata recording is done on a best-effort basis.  If metadata recording fails for any reason, a warning is logged, but no error is raised.

"""

import logging

_log = logging.getLogger(__name__)

