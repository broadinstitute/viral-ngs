'''This gives a main() function that serves as a nice wrapper
around other commands and presents the ability to serve up multiple
command-line functions from a single python script.
'''

import os, os.path, tempfile, sys, shutil, logging, argparse
import util.version

__author__ = "dpark@broadinstitute.org"
__version__ = util.version.get_version()

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
    return os.path.basename(sys.argv[0]).rsplit('.',1)[0]

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
                help="Base directory for temp files. [default: %(default)s]",
                default=v)
            parser.add_argument("--tmpDirKeep",
                action="store_true", dest="tmpDirKeep",
                help="""Keep the tmpDir if an exception occurs while
                    running. Default is to delete all temp files at
                    the end, even if there's a failure.""",
                default=False)
        elif k=='version':
            if not v:
                v=__version__
            parser.add_argument('--version', '-V', action='version', version=v)
        else:
            raise Exception("unrecognized argument %s" % arg)
    return parser

def main_command(mainfunc):
    ''' This wraps a python method in another method that can be called
        with an argparse.Namespace object. When called, it will pass all
        the values of the object on as parameters to the function call.
    '''
    def _main(args):
        args2 = dict((k,v) for k,v in vars(args).items() if k not in ('loglevel','tmpDir','tmpDirKeep','version','func_main','command'))
        mainfunc(**args2)
    _main.__doc__ = mainfunc.__doc__
    return _main

def attach_main(parser, cmd_main, split_args=False):
    ''' This attaches the main function call to a parser object.
    '''
    if split_args:
        cmd_main = main_command(cmd_main)
    parser.description = cmd_main.__doc__
    parser.set_defaults(func_main=cmd_main)
    return parser

def make_parser(commands, description):
    ''' commands: a list of pairs containing the following:
            1. name of command (string, no whitespace)
            2. method to call (no arguments) that returns an argparse parser.
            If commands contains exactly one member and the name of the
            only command is None, then we get rid of the whole multi-command
            thing and just present the options for that one function.
        description: a long string to present as a description of your script
            as a whole if the script is run with no arguments
    '''
    if len(commands)==1 and commands[0][0]==None:
        # only one (nameless) command in this script, simplify
        parser = commands[0][1]()
        parser.set_defaults(command='')
    else:
        # multiple commands available
        parser = argparse.ArgumentParser(description=description,
            usage='%(prog)s subcommand', add_help=False)
        parser.add_argument('--help', '-h', action='help', help=argparse.SUPPRESS)
        parser.add_argument('--version', '-V', action='version', version=__version__, help=argparse.SUPPRESS)
        subparsers = parser.add_subparsers(title='subcommands', dest='command')
        for cmd_name, cmd_parser in commands:
            p = subparsers.add_parser(cmd_name)
            cmd_parser(p)
    return parser

def main_argparse(commands, description):
    parser = make_parser(commands, description)
    
    # if called with no arguments, print help
    if len(sys.argv)==1:
        parser.parse_args(['--help'])
    elif len(sys.argv)==2 and (len(commands)>1 or commands[0][0]!=None):
        parser.parse_args([sys.argv[1], '--help'])
    args = parser.parse_args()
    
    setup_logger(not hasattr(args, 'loglevel') and 'DEBUG' or args.loglevel)
    log.info("software version: %s, python version: %s" % (__version__, sys.version))
    log.info("command: %s %s %s" % (
        sys.argv[0], sys.argv[1],
        ' '.join(["%s=%s" % (k,v) for k,v in vars(args).items() if k not in ('command', 'func_main')])))
    
    if hasattr(args, 'tmpDir'):
        ''' If this command has a tmpDir option, use that as a base directory
            and create a subdirectory within it which we will then destroy at
            the end of execution.
        '''
        proposed_dir = 'tmp-%s-%s' % (script_name(),args.command)
        if 'LSB_JOBID' in os.environ:
            proposed_dir = 'tmp-%s-%s-%s-%s' % (script_name(),args.command,os.environ['LSB_JOBID'],os.environ['LSB_JOBINDEX'])
        tempfile.tempdir = tempfile.mkdtemp(prefix='%s-'%proposed_dir, dir=args.tmpDir)
        log.debug("using tempDir: %s" % tempfile.tempdir)
        os.environ['TMPDIR'] = tempfile.tempdir     # this is for running R
        try:
            ret = args.func_main(args)
        except:
            if hasattr(args, 'tmpDirKeep') and args.tmpDirKeep and not (tempfile.tempdir.startswith('/tmp') or tempfile.tempdir.startswith('/local')):
                log.exception("Exception occurred while running %s, saving tmpDir at %s" % (args.command, tempfile.tempdir))
            else:
                shutil.rmtree(tempfile.tempdir)
            raise
        else:
            shutil.rmtree(tempfile.tempdir)
    else:
        # otherwise just run the command
        ret = args.func_main(args)
    if ret==None:
        ret = 0
    return ret

def find_tmpDir():
    ''' This provides a suggested base directory for a temp dir for use in your
        argparse-based tmpDir option.
    '''
    tmpdir = '/tmp'
    if os.access('/local/scratch', os.X_OK | os.W_OK | os.R_OK) and os.path.isdir('/local/scratch'):
        tmpdir = '/local/scratch'
    if 'LSB_JOBID' in os.environ:
        # this directory often exists for LSF jobs, but not always.
        # for example, if the job is part of a job array, this directory is called
        # something unpredictable and unfindable, so just use /local/scratch
        proposed_dir = '/local/scratch/%s.tmpdir' % os.environ['LSB_JOBID']
        if os.access(proposed_dir, os.X_OK | os.W_OK | os.R_OK):
            tmpdir = proposed_dir
    elif 'TMPDIR' in os.environ and os.path.isdir(os.environ['TMPDIR']):
        tmpdir = os.environ['TMPDIR']
    return tmpdir
