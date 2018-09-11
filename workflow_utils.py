#!/usr/bin/env python
"""
Utilities for dealing with workflows, including cloud workflows.

Commands to help jointly version code and data using tools such as git-annex and DataLad.
"""

__author__ = "ilya@broadinstitute.org"
__commands__ = []

import platform
assert platform.python_version().startswith('2.7')  # dxpy requirement; should drop and use command-line

import argparse
import logging
import json
import subprocess
import os
import os.path
import shutil
import glob
import collections
import time
import getpass
import uuid
import SimpleHTTPServer
import SocketServer

import dxpy
import dxpy.bindings.dxfile_functions
import boto3

import util.cmd
import util.file
import util.misc

_log = logging.getLogger(__name__)

logging.basicConfig(format="%(asctime)s - %(module)s:%(lineno)d:%(funcName)s - %(levelname)s - %(message)s")
_log.setLevel(logging.DEBUG)

def _is_str(obj):
    """Test if obj is a string type, in a python2/3 compatible way.
    From https://stackoverflow.com/questions/4232111/stringtype-and-nonetype-in-python3-x
    """
    try:
        return isinstance(obj, basestring)
    except NameError:
        return isinstance(obj, str)

def _pretty_print_json(json_dict):
    """Return a pretty-printed version of a dict converted to json, as a string."""
    return json.dumps(json_dict, indent=4, separators=(',', ': '))

def _run(cmd):
    print('running command: ', cmd)
    beg_time = time.time()
    subprocess.check_call(cmd, shell=True)
    print('command succeeded in {}s: {}'.format(time.time()-beg_time, cmd))

def _run_get_output(cmd):
    print('running command: ', cmd)
    beg_time = time.time()
    output = subprocess.check_output(cmd, shell=True)
    print('command succeeded in {}s: {}'.format(time.time()-beg_time, cmd))
    return output

def _run_get_json(cmd):
    return json.loads(_run_get_output(cmd).strip())

def workflow_utils_init():
    """Install the dependencies: cromwell and dxpy and git-annex."""
    _run('conda install cromwell dxpy git-annex')

#
# given a dx analysis, run it locally
#

def get_workflow_inputs(workflow_name, t_dir, docker_img):
    """Run womtool to get the inputs of the wdl workflow"""
    with util.file.pushd_popd(t_dir):
        _run('docker run --rm ' + docker_img + ' tar cf - source/pipes/WDL > wdl.tar')
        _run('tar xvf wdl.tar')
        for f in glob.glob('source/pipes/WDL/workflows/*.wdl'):
            shutil.copy(f, '.')
        for f in glob.glob('source/pipes/WDL/workflows/tasks/*.wdl'):
            shutil.copy(f, '.')
        shutil.rmtree('source')
        os.unlink('wdl.tar')
        return _run_get_json('womtool inputs ' + workflow_name + '.wdl')
    

def _get_dx_val(val, dx_files_dir = 'input_files'):
    """Resolve a dx value: if it is a scalar, just return that;
    if it is a dx file, fetch the file, cache it locally, and return the path to the file."""

    print('parsing val: ', val)
    util.file.mkdir_p(dx_files_dir)
    if isinstance(val, collections.Mapping) and '$dnanexus_link' in val:
        link = val['$dnanexus_link']
        if isinstance(link, collections.Mapping):
            return _get_dx_val(dxpy.DXAnalysis(link['analysis']).describe()['output'][link['stage']+'.'+link['field']])
        elif link.startswith('file-'):
            dxid = link
            descr = dxpy.describe(dxid)
            dx_file = os.path.join(dx_files_dir, dxid) + '-' + descr['name']
            file_size = int(descr['size'])

            # see if the file is cached in git-annex
            ga_mdata = _run_get_json('git annex metadata --json --key=WORM-s0-m0--dx-' + dxid)
            if ga_mdata['fields']:
                ga_key = ga_mdata['fields']['ga_key'][0]
                _run('git annex get --key ' + ga_key)
                _run('git annex fromkey ' + ga_key + ' ' + dx_file)

            if not os.path.isfile(dx_file) or os.path.getsize(dx_file) != file_size:
                print('fetching', dxid, 'to', dx_file)
                # TODO: check that there is enough free space (with some to spare)
                if os.path.isfile(dx_file):
                    os.unlink(dx_file)
                fs_info = os.statvfs(dx_files_dir)
                assert file_size < fs_info.f_bsize * fs_info.f_bavail
                dxpy.bindings.dxfile_functions.download_dxfile(dxid=dxid, filename=dx_file+'.fetching', show_progress=True)
                assert os.path.getsize(dx_file+'.fetching') == file_size
                os.rename(dx_file+'.fetching', dx_file)
                print('fetched', dxid, 'to', dx_file)
                print('curdir is ', os.getcwd())
                _run('git annex add ' + dx_file)
                ga_key = _run_get_output('git annex lookupkey ' + dx_file).strip()
                _run('git annex metadata -s dxid+=' + dxid + ' ' + dx_file)
                # record a mapping from the dxid to the git-annex key
                _run('git annex metadata --key=WORM-s0-m0--{} -s ga_key={}'.format('dx-'+dxid, ga_key))
                # register a URL that can be used to re-fetch this file from DNAnexus;
                # the URL is only valid if the 'run_dx_url_server' command is running.
                _run('git annex registerurl ' + ga_key + ' ' + ' http://localhost:8080/dx/' + dxid)
        else:
            raise RuntimeError('Unknown dxid {}'.format(dxid))
        return dx_file
    else:
        return val
# end: def _get_dx_val(val, dx_files_dir = 'input_files')

def _parse_cromwell_output_str(cromwell_output_str):
    """Parse cromwell output"""
    assert cromwell_output_str.count('Final Outputs:') == 1
    json_beg = cromwell_output_str.index('Final Outputs:') + len('Final Outputs:')
    json_end = cromwell_output_str.index('\n}\n', json_beg) + 2
    return json.loads(cromwell_output_str[json_beg:json_end])

def run_analysis_wdl(workflow_name, analysis_inputs_from_dx_analysis, docker_img, analysis_dir, analysis_inputs_specified=None, 
                     analysis_descr='', analysis_labels=None,
                     data_repo=None, data_remote=None):
    """Run a WDL analysis locally, using a specified version of viral-ngs.

    Notes:

       - record the actual git hash (since git tag might change)

    """
    assert not os.path.isabs(analysis_dir)

    if not analysis_inputs_specified:
        analysis_inputs_specified = {}
    else:
        print('analysis_inputs_specified=', util.file.slurp_file(analysis_inputs_specified).strip())
        analysis_inputs_specified = json.loads(util.file.slurp_file(analysis_inputs_specified).strip())

    analysis_id = create_analysis_id(workflow_name)
    print('ANALYSIS_ID is ', analysis_id)

    _run('docker pull ' + docker_img)
    docker_img_hash = docker_img + '@' + get_docker_hash(docker_img)

    dx_analysis = dxpy.DXAnalysis(dxid=analysis_inputs_from_dx_analysis)
    dx_analysis_descr = dx_analysis.describe()

    assert os.path.exists('.git/annex')
    data_repo = data_repo or os.getcwd()

    with util.file.tmp_dir('_workflow_copy') as t_dir_orig:
        with util.file.pushd_popd(t_dir_orig):
            _run('git clone ' + data_repo + ' data')
            t_dir_git = os.path.join(t_dir_orig, 'data')
            util.file.mkdir_p(t_dir_git)
            with util.file.pushd_popd(t_dir_git):
                _run('git annex init')
                _run('git annex describe here running_' + analysis_id)
                _run('git config annex.security.allowed-http-addresses "127.0.0.1 ::1 localhost"')
                if data_remote:
                    _run('git annex enableremote ' + data_remote)

                t_dir = os.path.join(t_dir_git, analysis_dir + '-' + analysis_id)
                util.file.mkdir_p(t_dir)
                with util.file.pushd_popd(t_dir):
                    print('TTTTTTTTTTT t_dir=', t_dir)

                    output_dir = 'output'

                    util.file.mkdir_p('input_files')

                    wdl_wf_inputs = get_workflow_inputs(workflow_name, t_dir, docker_img=docker_img_hash)

                    # TODO: use git annex batch mode to determine the keys for all the file-xxxx files, then
                    # use batch mode to get the keys.  use -J to parallelize.  also option to use a cache remote, as
                    # described at https://git-annex.branchable.com/tips/local_caching_of_annexed_files/

                    # TODO: diff stages may have same-name inputs
                    dx_wf_inputs = { k.split('.')[-1] : v for k, v in dx_analysis_descr['runInput'].items()}   # check for dups
                    dx_wf_orig_inputs = { k.split('.')[-1] : v for k, v in dx_analysis_descr['originalInput'].items()}   # check for dups
                    print('DX_WF_RUNINPUTS', '\n'.join(dx_wf_inputs.keys()))
                    print('DX_WF_ORIGINPUTS', '\n'.join(dx_wf_orig_inputs.keys()))
                    new_wdl_wf_inputs = {}
                    for wdl_wf_input, wdl_wf_input_descr in wdl_wf_inputs.items():
                        wdl_wf_input_full = wdl_wf_input
                        wdl_wf_input = wdl_wf_input.split('.')[-1]
                        if wdl_wf_input in analysis_inputs_specified:
                            dx_wf_input = analysis_inputs_specified[wdl_wf_input]
                            new_wdl_wf_inputs[wdl_wf_input_full] = map(_get_dx_val, dx_wf_input) \
                                                                   if isinstance(dx_wf_input, list) \
                                                                      else _get_dx_val(dx_wf_input)
                        elif wdl_wf_input in dx_wf_inputs:
                            print('HAVE', wdl_wf_input, wdl_wf_input_descr, dx_wf_inputs[wdl_wf_input])
                            dx_wf_input = dx_wf_inputs[wdl_wf_input]
                            new_wdl_wf_inputs[wdl_wf_input_full] = map(_get_dx_val, dx_wf_input) \
                                                                   if isinstance(dx_wf_input, list) \
                                                                      else _get_dx_val(dx_wf_input)
                        elif '(optional' not in wdl_wf_input_descr and wdl_wf_input in dx_wf_orig_inputs:
                            print('HAVE', wdl_wf_input, wdl_wf_input_descr, dx_wf_orig_inputs[wdl_wf_input])
                            dx_wf_input = dx_wf_orig_inputs[wdl_wf_input]
                            new_wdl_wf_inputs[wdl_wf_input_full] = map(_get_dx_val, dx_wf_input) \
                                                                   if isinstance(dx_wf_input, list) \
                                                                      else _get_dx_val(dx_wf_input)
                        else:
                            print('MISSING', wdl_wf_input, wdl_wf_input_descr)
                            assert '(optional' in wdl_wf_input_descr


                    print(_pretty_print_json(new_wdl_wf_inputs))

                    ################# put in the right docker ID!!  and find a place to keep the docker cache.

                    with open('inputs.json', 'wt') as wf_out:
                        json.dump(new_wdl_wf_inputs, wf_out, indent=4, separators=(',', ': '))

                    # TODO: option to update just some of the tasks.
                    # actually, when compiling WDL, should have this option -- or, actually,
                    # should make a new workflow where older apps are reused for stages that have not changed.
                    _run('sed -i -- "s|{}|{}|g" *.wdl'.format('quay.io/broadinstitute/viral-ngs', docker_img_hash))

                    util.file.mkdir_p(output_dir)
                    util.file.mkdir_p(os.path.join(output_dir, 'outputs'))
                    util.file.mkdir_p(os.path.join(output_dir, 'logs'))
                    util.file.mkdir_p(os.path.join(output_dir, 'call_logs'))
                    #util.file.mkdir_p(os.path.join(output_dir, 'metadata'))
                    wf_opts_dict = { "final_workflow_outputs_dir": os.path.join(output_dir, 'outputs'),
                                     "final_workflow_log_dir": os.path.join(output_dir, 'logs'),
                                     "final_call_logs_dir": os.path.join(output_dir, 'call_logs')
                    }
                    util.file.dump_file('cromwell_opts.json', _pretty_print_json(wf_opts_dict))
                    util.file.dump_file('execution_env.json', _pretty_print_json(dict(ncpus=util.misc.available_cpu_count())))

                    util.file.dump_file('analysis_labels.json',
                                        _pretty_print_json(dict(analysis_descr=analysis_descr,
                                                                docker_img=docker_img,
                                                                docker_img_hash=docker_img_hash,
                                                                inputs_from_dx_analysis=analysis_inputs_from_dx_analysis,
                                                                **dict(analysis_labels or {}))))

                    # add cromwell labels: dx project, the docker tag we ran on, etc.

                    _log.info('Validating workflow')
                    _run('womtool validate -i inputs.json ' + workflow_name + '.wdl')
                    _log.info('Validated workflow; calling cromwell')
                    try:
                        cromwell_output_str = _run_get_output('cromwell run ' + workflow_name + \
                                                              '.wdl -i inputs.json -l analysis_labels.json ' + \
                                                              ' -o cromwell_opts.json' + \
                                                              ' -m ' + \
                                                              os.path.join(output_dir, 'cromwell_execution_metadata.json'))
                        cromwell_returncode = 0
                    except subprocess.CalledProcessError as called_process_error:
                        cromwell_output_str = called_process_error.output
                        cromwell_returncode = called_process_error.returncode
                        
                    _log.info('Cromwell returned with return code %d', cromwell_returncode)

                    util.file.dump_file(os.path.join(output_dir, 'cromwell_output.txt'), cromwell_output_str)

                    if cromwell_returncode == 0:
                        def make_paths_relative(v):
                            print('make_paths_relative: v=', v, 't_dir_git=', t_dir_git)
                            if _is_str(v) and os.path.isabs(v) and v.startswith(t_dir):
                                return os.path.relpath(v, t_dir)
                            if isinstance(v, list):
                                return list(map(make_paths_relative, v))
                            return v
                        cromwell_output_json = {k: make_paths_relative(v)
                                                for k, v in _parse_cromwell_output_str(cromwell_output_str).items()}
                        util.file.dump_file('outputs.json',
                                            _pretty_print_json(cromwell_output_json))
                        util.file.make_empty('analysis_succeeded.txt')
                    else:
                        util.file.make_empty('analysis_failed.txt')

                    _run('rm *.wdl')

                _run('sudo chown -R $USER . || true')
                _run('git annex add')
                _run('git commit -m after_running_analysis_' + analysis_id + '.')
#                _run('git annex initremote content type=directory directory=/ndata/git-annex-content/ encryption=none')
                _run('git annex move --all --to {} -J{}'.format(data_remote or 'origin', util.misc.available_cpu_count()))
                _run('git annex dead here')
                _run('git annex sync --message git_annex_sync_analysis_{}'.format(analysis_id))

                # enable cleanup
                _run('chmod -R u+w . || true')


def create_analysis_id(workflow_name):
    """Generate a unique ID for the analysis."""
    return util.file.string_to_file_name('-'.join(map(str, 
                                                       ('analysis', time.strftime('%Y%m%d-%H%M%S', time.localtime())[2:], 
                                                        uuid.uuid4(), workflow_name))))[:1024]

def get_docker_hash(docker_img):
    if docker_img.startswith('sha256:'):
        return docker_img
    digest_line = _run_get_output('docker images ' + docker_img + ' --digests --no-trunc --format '
                                  '"{{.Repository}}:{{.Tag}} {{.Digest}}"')
    assert digest_line.count('\n') == 1
    img, digest = digest_line.strip().split()
    assert img == docker_img + (':latest' if ':' not in docker_img else '')
    assert digest.startswith('sha256:') and len(digest) == 71
    return digest

def hash_and_store_file(fname):
    '''Store the given file in a key-value store, '''
    pass

def gather_run_results(docker_img, cromwell_output):
    pass
    
    #asdf

    # so, we also need to generate a unique id for this run;
    # since there can be multiple runs.

    # also, need to see how to handle failures, full or partial.
    #
#     sdb_client = boto3.client('sdb')
#    response = sdb_client.put_attributes(DomainName='viral_ngs_benchmarks',
#                                         ItemName=create_run_id(),
#                                         Attributes=[
                                            
#                                         ])

# =======================

########################################################################################################################


def parser_run_analysis_wdl(parser=argparse.ArgumentParser()):
    parser.add_argument('workflow_name', help='Workflow name')
    parser.add_argument('analysis_dir', default='runs/an',
                        help='directory where analysis will be stored; a unique suffix will be added')
    parser.add_argument('--analysisInputsFromDxAnalysis', dest='analysis_inputs_from_dx_analysis',
                        help='DNAnexus analysis ID to take analysis inputs from; specific ones can be overridden by '
                        '--analysisInputsSpecified')
    parser.add_argument('--dockerImg', dest='docker_img', default='quay.io/broadinstitute/viral-ngs')
    parser.add_argument('--analysisInputsSpecified', dest='analysis_inputs_specified',
                        help='explicitly specified analysis inputs')
    parser.add_argument('--analysisDescr', dest='analysis_descr', help='description of the run')
    parser.add_argument('--dataRepo', dest='data_repo', help='git data repository')
    parser.add_argument('--dataRemote', dest='data_remote', help='git-annex data remote')
    parser.add_argument('--analysisLabels', dest='analysis_labels', nargs=2, action='append',
                        help='labels to attach to the analysis')

    util.cmd.attach_main(parser, run_analysis_wdl, split_args=True)

__commands__.append(('run_analysis_wdl', parser_run_analysis_wdl))


########################################################################################################################

class RedirectToDNAnexus(SimpleHTTPServer.SimpleHTTPRequestHandler):

   def do_GET(self):
       try:
           print('path=', self.path)
           try:
               file_info = _run_get_json('dx describe --json ' + self._get_dx_id())
           except subprocess.CalledProcessError as e:
               print('CalledProcessError:', e)
               raise
           print('got file_info', file_info)
           assert file_info['id'] == self._get_dx_id()

           self.send_response(307)
           if False and 'media' in file_info:
               self.send_header('Content-Type', file_info['media'])
           if False and 'size' in file_info:
               self.send_header('Content-Length', str(file_info['size']))
           #self.send_header('Content-Disposition', 'attachment; filename="{}"'.format(file_info['name']))

           new_path = _run_get_output('dx make_download_url --duration 2h ' + self._get_dx_id())
           print('new_path=', new_path)
           self.send_header('Location', new_path)
           self.end_headers()
       except (IOError, subprocess.CalledProcessError):
           self.send_error(404, 'file not found')

   def do_HEAD(self):
       self.do_GET()

   def _get_dx_id(self):
       return self.path[len('/dx/'):]

def run_dx_url_server(dummy, port=8080):
    """Start a webserver that will redirect URL requests of the form http://localhost/dx/file-xxxxxx to DNAnexus.
    This gives each dx file a stable URL, permitting such files to be added to git-annex with the addurl command."""

    server = None
    try:           
        SocketServer.TCPServer.allow_reuse_address = True
        server = SocketServer.TCPServer(("", port), RedirectToDNAnexus)
        server.serve_forever()
    except Exception as e:
        if server is None:
            print('No server!')
        else:
            print('calling shutdown...')
            server.shutdown()
            print('calling server_close...')
            server.server_close()
            print('re-raising exception', e)
            raise

def parser_run_dx_url_server(parser=argparse.ArgumentParser()):
    parser.add_argument('dummy', help='Ignored argument (running command with no args just prints help)')
    parser.add_argument('--port', default=8080, help='Port on which to run the webserver')
    util.cmd.attach_main(parser, run_dx_url_server, split_args=True)

__commands__.append(('run_dx_url_server', parser_run_dx_url_server))

########################################################################################################################

def full_parser():
    return util.cmd.make_parser(__commands__, __doc__)


if __name__ == '__main__':
    #print(_pretty_print_json(_parse_cromwell_output_str(util.file.slurp_file('/dev/shm/cromwell_out/wdl_output.txt'))))
    #print(get_docker_hash('quay.io/broadinstitute/viral-ngs'))
    #record_run_results(0,0)
    if False:
        run_dx_locally(workflow_name='assemble_denovo', analysis_dxid='analysis-FJfqjg005Z3Vp5Q68jxzx5q1',
                       docker_img='quay.io/broadinstitute/viral-ngs')
    util.cmd.main_argparse(__commands__, __doc__)
