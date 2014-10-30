'''This gives a main() function that serves as a nice wrapper
around other commands and presents the ability to serve up multiple
command-line functions from a single python script.
'''

__author__ = "dpark@broadinstitute.org"
__version__ = "PLACEHOLDER"
__date__ = "PLACEHOLDER"

import os, tempfile, sys, shutil, logging

log = logging.getLogger()
tmpDir = None

def setup_logger(log_level):
    loglevel = getattr(logging, log_level.upper(), None)
    assert loglevel, "unrecognized log level: %s" % log_level
    log.setLevel(loglevel)
    h = logging.StreamHandler()
    h.setFormatter(logging.Formatter("%(asctime)s - %(module)s:%(lineno)d:%(funcName)s - %(levelname)s - %(message)s"))
    log.addHandler(h)

def script_name():
    return sys.argv[0].split('/')[-1].rsplit('.',1)[0]

def common_args(parser, arglist=(('tmpDir',None), ('loglevel',None))):
    for k,v in arglist:
        if k=='loglevel':
            if not v:
                v = 'DEBUG'
            parser.add_argument("--loglevel", dest="loglevel",
                help="Verboseness of output.  [default: %(default)s]",
                default=v,
                choices=('DEBUG','INFO','WARNING','ERROR','CRITICAL','EXCEPTION'))
        elif k=='tmpDir':
            if not v:
                v = find_tmpDir()
            parser.add_argument("--tmpDir", dest="tmpDir",
                help="Directory for temp files.  [default: %(default)s]",
                default=v)
        elif k=='tmpDirKeep':
            if v==None:
                v = False
            else:
                v = v and True or False
            parser.add_argument("--tmpDirKeep",
                action="store_true", dest="tmpDirKeep",
                help="Delete the tmpDir, even if an exception occurs while running. [default: %(default)s]",
                default=v)
            parser.add_argument("--noTmpDirKeep",
                action="store_false", dest="tmpDirKeep",
                help="Delete the tmpDir if an exception occurs while running. [default: %(default)s]")
        elif k=='version':
            if not v:
                v=__version__
            parser.add_argument("--version", action='version', version=v)
        else:
            raise Exception("unrecognized argument %s" % arg)
    return parser

def main_argparse(commands, description):
    ''' commands: a list of 3-tuples containing the following:
            1. name of command (string, no whitespace)
            2. method to call that takes one argument (an argparse construct),
                and returns the desired exit code
            3. method to call (no arguments) that returns an argparse parser.
            If commands contains exactly one member and the name of the
            only command is None, then we get rid of the whole multi-command
            thing and just present the options for that one function.
        tool_paths: a dict.  we will set the 'tmpDir' value so that your
            commands will have access to a suggested temp directory
        description: a long string to present as a description of your script
            as a whole if the script is run with no arguments
    '''
    assert description, "docstring cannot be absent!"
    tmpDir = find_tmpDir()

    cmdlist = [x[0] for x in commands]
    commands = dict([(x[0],x[1:]) for x in commands])

    if len(cmdlist)==1 and cmdlist[0]==None:
        # only one (nameless) command in this script, simplify
        command = None
        argv = sys.argv[1:]
    else:
        # multiple commands available
        if len(sys.argv) <= 1:
            print ("Usage: python %s commandname options" % sys.argv[0])
            if description.strip():
                print (description)
            print ("\ncommands:")
            for cmd in cmdlist:
                print ("\t%s" % cmd)
            print ("\nRun a command with no options for help on that command.")
            return
        command = sys.argv[1]
        assert command in commands, "command '%s' not recognized" % command
        argv = sys.argv[2:]

    parser = commands[command][1]()
    if len(argv)==0:
        parser.print_help()
        return
    args = parser.parse_args(argv)

    setup_logger(not hasattr(args, 'loglevel') and 'DEBUG' or args.loglevel)
    log.info("software version: " + __version__)
    log.debug("python version: " + sys.version)
    log.debug("command line parameters (including implicit defaults): %s" % (
        ' '.join(["%s=%s" % (k,v) for k,v in vars(args).items()])))

    if hasattr(args, 'tmpDir'):
        ''' If this command has a tmpDir option, use that as a base directory
            and create a subdirectory within it which we will then destroy at
            the end of execution.
        '''
        proposed_dir = 'tmp-%s-%s' % (script_name(),command!=None and command or '')
        if 'LSB_JOBID' in os.environ:
            proposed_dir = 'tmp-%s-%s-%s-%s' % (script_name(),command,os.environ['LSB_JOBID'],os.environ['LSB_JOBINDEX'])
        tempfile.tempdir = tempfile.mkdtemp(prefix='%s-'%proposed_dir, dir=args.tmpDir)
        log.debug("using tempDir: %s" % tempfile.tempdir)
        os.environ['TMPDIR'] = tempfile.tempdir     # this is for running R
        try:
            ret = commands[command][0](args)
        except:
            if hasattr(args, 'tmpDirKeep') and args.tmpDirKeep and not (tempfile.tempdir.startswith('/tmp') or tempfile.tempdir.startswith('/local')):
                log.exception("Exception occurred while running %s, saving tmpDir at %s" % (command, tempfile.tempdir))
            else:
                shutil.rmtree(tempfile.tempdir)
            raise
        else:
            shutil.rmtree(tempfile.tempdir)
        return ret
    else:
        # otherwise just run the command
        return commands[command][0](args)

def find_tmpDir():
    ''' This provides a suggested base directory for a temp dir for use in your
        optparse-based tmpDir option.
    '''
    tmpdir = '/tmp'
    if os.access('/local/scratch', os.X_OK | os.W_OK | os.R_OK):
        tmpdir = '/local/scratch'
    if 'LSB_JOBID' in os.environ:
        # this directory often exists for LSF jobs, but not always.
        # for example, if the job is part of a job array, this directory is called
        # something unpredictable and unfindable, so just use /local/scratch
        proposed_dir = '/local/scratch/%s.tmpdir' % os.environ['LSB_JOBID']
        if os.access(proposed_dir, os.X_OK | os.W_OK | os.R_OK):
            tmpdir = proposed_dir
    return tmpdir
