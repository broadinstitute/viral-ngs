import collections
import os
import os.path
import json
import time

from util._metadata.file_arg import FileArg

import networkx
import networkx.algorithms.dag

# ** class ProvenanceGraph
class ProvenanceGraph(networkx.DiGraph):

    '''Provenance graph representing a set of computations.  It has two types of nodes: steps and files.  A step node (snode)
    represents a computation step (a specific execution of a specific command); a file node (fnode) represents a particular file
    that was read/written by one or more steps.

    Edges go from input files of a step to the step, and from a step to its output files, so the graph is bipartite.
    Various information about the parameters and execution of a step is represented as attributes of the step node.
    '''

    FileNode = collections.namedtuple('FileNode', 'realpath file_hash mtime')

    def __init__(self):
        """Initialize an empty provenance graph."""
        super(ProvenanceGraph, self).__init__()

    def is_fnode(self, n): return self.nodes[n]['node_kind'] == 'file'
    def is_snode(self, n): return self.nodes[n]['node_kind'] == 'step'

    @property
    def file_nodes(self): 
        """The nodes representing files"""
        return [n for n in self if self.is_fnode(n)]


    @property
    def step_nodes(self): 
        """The nodes representing steps"""
        return [n for n in self if self.is_snode(n)]

# *** Loading pgraph
    def load(self, metadata_path=None, max_age_days=10000):
        """Read the provenance graph from the metadata directory.

        Args:
           path: path to the metadata directory; if None (default), use the path specified by the environment var VIRAL_NGS_METADATA_PATH.
           max_age_days: ignore steps older than this many days
        """

        #
        # This is not quite a straightforward "load" operation -- we transform the data a bit in the process.
        # When recording, we record everything we might later need; but when loading, we might filter out some info.
        # We might also compute some useful derived info.
        #

        metadata_path = metadata_path or metadata_dir()

# **** Load steps and files
        for step_record_fname in os.listdir(metadata_path or metadata_dir()):
            if not step_record_fname.endswith('.json'): continue

            _log.info('loading step {}'.format(step_record_fname))
            json_str = util.file.slurp_file(os.path.join(metadata_path, step_record_fname))

            step_record = json.loads(json_str)

            # check that all expected elements are in the json file -- some older-format files might be missing data
            if not is_valid_step_record(step_record): 
                _log.warning('not valid step: {}'.format(step_record_fname))
                continue
            if step_record['step']['run_info']['exception']: continue  # for now, ignore steps that failed
            if step_record['step'].get('enclosing_steps', ''): continue  # for now, skip steps that are sub-steps of other steps
            if ((time.time() - step_record['step']['run_info']['beg_time']) / (60*60*24)) > max_age_days: continue

            step_id = step_record['step']['step_id']

            # fix an issue in some legacy files
            step_record['step']['metadata_from_cmd_line'] = dict(step_record['step'].get('metadata_from_cmd_line', {})) 
            step_record['step']['step_name'] = step_record['step']['metadata_from_cmd_line'].get('step_name',
                                                                                                 step_record['step']['cmd_name'])

            self.add_node(step_id, node_kind='step', step_record_fname=step_record_fname, **step_record['step'])

            # Add input and output files of this step as data nodes.
            for arg, val in step_record['step']['args'].items():
                def gather_files(val):
                    """Return a flattened list of files denoted by this command argument (or empty list if it does not denote files)."""
                    return (FileArg.is_from_dict(val) and [val]) \
                        or (isinstance(val, (list, tuple)) and functools.reduce(operator.concat, list(map(gather_files, val)), [])) or []

                for files in gather_files(val):
                    assert FileArg.is_from_dict(files)
                    for f in files['files']:
                        if dict_has_keys(f, 'hash fname realpath size'):
                            file_node = self.FileNode(f['realpath'], f['hash'], f['mtime'])
                            self.add_node(file_node, node_kind='file', **f)
                            e = (file_node, step_id) if files['mode'] == 'r' else (step_id, file_node)
                            self.add_edge(*e, arg=arg)

                            # gather any per-file metadata specified on the command line.  currently, such metadata can be specified
                            # for a given command arg, and applied to all files denoted by the arg.
                            for metadata_attr, metadata_val in step_record['step']['metadata_from_cmd_line'].items():
                                pfx = 'file.{}.'.format(arg)
                                if metadata_attr.startswith(pfx):
                                    self[e[0]][e[1]][metadata_attr[len(pfx):]] = metadata_val
                    # end: for f in files['files']
                # end: for files in gather_files(val)
            # end: for arg, val in step_record['step']['args'].items()
        # end: for step_record_fname in os.listdir(path)

# **** Check for anomalies

        for f in self.file_nodes:
            assert self.in_degree[f] <= 1
            # if self.in_degree[f] > 1:
            #     _log.warning('ANOMALY: file with indegree {}: {}'.format(self.in_degree[f], f))
            #     for e in self.in_edges(nbunch=f, data=True):
            #         _log.warning(e[0], e[1], json.dumps(e[2], indent=4))

        self._reconstruct_missing_connections()

        assert networkx.algorithms.dag.is_directed_acyclic_graph(self)

        #
        # Print nodes with unknown origin
        # 

        unknown_origin_files = []
        for file_node in self.file_nodes:
            #assert self.in_degree[file_node] <= 1
            if self.in_degree[file_node] > 1:
                _log.warning('indegree of {} is {}'.format(file_node, self.in_degree[file_node]))
            if self.in_degree[file_node] == 0 and os.path.isfile(file_node.realpath):
                unknown_origin_files.append(file_node.realpath)
                #_log.warning('UNKNOWN ORIGIN: {}'.format(file_node))

        _log.warning('UNKNOWN_ORIGIN:\n{}'.format('\n'.join(unknown_origin_files)))

    # end: def load(self, path=None)

    def _reconstruct_missing_connections(self):
        """Reconstruct missing connections between steps and files.
        """

        hash2files, realpath2files = self._compute_hash_and_realpath_indices()

        for f in list(self.file_nodes):
            if self.is_fnode(f) and not self.pred[f]:
                # file f is read by some step, but we don't have a record of a step that wrote this exact file.
                # do we have 
                print('trying to reconnect to {}'.format(f))
                f2s = [ f2 for f2 in (hash2files[f.file_hash] & realpath2files[f.realpath]) if f2 != f and f2.mtime < f.mtime ]
                print('f2s={}'.format(f2s))
                if f2s:
                    f2s = sorted(f2s, key=operator.attrgetter('mtime'))
                    for s in list(self.succ[f]):
                        f2s_bef = [ f2 for f2 in f2s if f2.mtime < self.nodes[s]['run_info']['beg_time'] ]
                        if f2s_bef:
                            f2 = f2s_bef[-1]
                            assert f2 != f and f2.file_hash == f.file_hash and f2.mtime < f.mtime
                            self.add_edge(f2, s, **self.edges[f, s])
                            self.remove_edge(f, s)
                            
                            _log.info('reconnected: {}->{}'.format(f2.realpath, s))
    # end: def _reconstruct_missing_connections(self, file_nodes=None):

    def _reconstruct_fnode_maker(self, fnode):
        """Reconstruct missing maker step of fnode: find another node (if exists) that represents the same file (in the same location)
        with the same contents, just with an earlier mtime.  Return the new fnode if found, else return none.
        """

        hash2files, realpath2files = self._compute_hash_and_realpath_indices()

        # file f is read by some step, but we don't have a record of a step that wrote this exact file.
        # do we have 
        f = fnode
        print('trying to reconnect to {}'.format(f))
        f2s = [ f2 for f2 in (hash2files[f.file_hash] & realpath2files[f.realpath]) if f2 != f ]
        print('f2s={}'.format(f2s))
        if f2s:
            f2s = sorted(f2s, key=operator.attrgetter('mtime'))
            f2 = f2s[-1]
            assert f2 != f and f2.file_hash == f.file_hash
            return f2
        return None
    # end: def _reconstruct_fnode_maker(self, fnode):

    def _compute_hash_and_realpath_indices(self):
        """Compute mapping from hash to file nodes and from realpath to file nodes"""

        hash2files = collections.defaultdict(set)
        realpath2files = collections.defaultdict(set)
        for f in self.file_nodes:
            hash2files[f.file_hash].add(f)
            realpath2files[f.realpath].add(f)

        for h, fs in hash2files.items():
            if len(fs) > 1:
                for f1 in fs:
                    for f2 in fs:
                        if f1.realpath != f2.realpath and os.path.isfile(f1.realpath) and os.path.isfile(f2.realpath) and \
                           os.path.samefile(f1.realpath, f2.realpath):
                            realpath2files[f1.realpath].add(f2)
                            realpath2files[f2.realpath].add(f1)
                            #_log.info('SAMEFILES:\n{}\n{}\n'.format(f1, f2))

        return hash2files, realpath2files
    # end: def _compute_hash_and_realpath_indices(self):

# *** write_dot
    def write_dot(self, dotfile, nodes=None, ignore_cmds=(), ignore_exts=(), title=''):
        """Write out this graph, or part of it, as a GraphViz .dot file."""

        ignore_exts = ()

        def get_val(d, keys):
            """Fetch a value from a nested dict using a sequence of keys"""
            if not keys: return d
            keys = util.misc.make_seq(keys)
            if isinstance(d, collections.Mapping) and keys[0] in d:
                return get_val(d[keys[0]], keys[1:])
            return None

        def fix_name(s, node2id={}):
            """Return a string based on `s` but valid as a GraphViz node name"""
            return 'n' + str(node2id.setdefault(s, len(node2id)))

        ignored = set()

        with open(dotfile, 'wt') as out:
            out.write('digraph G {\n')
            for n in self:
                if not nodes or n in nodes:
                    n_attrs = self.nodes[n]
                    if n_attrs['node_kind'] == 'step':
                        label = n_attrs['step_name']
                        if label in ignore_cmds: 
                            ignored.add(n)
                            continue
                        shape = 'invhouse'
                    else:
                        label = get_val(n_attrs, 'fname') or 'noname'
                        if any(label.endswith(e) for e in ignore_exts): 
                            ignored.add(n)
                            continue
                        shape = 'oval'
                        
                    out.write('{} [label="{}", shape={}];\n'.format(fix_name(n), label, shape))

            for u, v, arg in self.edges(data='arg'):
                if nodes and (u not in nodes or v not in nodes): continue
                if u in ignored or v in ignored: continue
                out.write('{} -> {} [label="{}"];\n'.format(fix_name(u), fix_name(v), arg))
            out.write('labelloc="t";\n')
            out.write('label="{}\n{}";\n'.format(time.strftime('%c'), title))
            out.write('}\n')
    # end: def write_dot(self, dotfile, nodes=None, ignore_cmds=(), ignore_exts=()):

    def write_svg(self, svgfile, *args, **kwargs):
        """Write out this graph to an .svg file.  See write_dot() for documentation of args."""
        with util.file.tempfname('.dot') as dotfile:
            self.write_dot(dotfile, *args, **kwargs)
            _shell_cmd('dot -Tsvg -o {} {}'.format(svgfile, dotfile), check=True)

# *** print_provenance
    def print_provenance(self, fname, svgfile=None):
        """Print the provenance info for the given file.
        
        Args:
            svgfile: if given, write the provenance graph for the given file to this .svg file
        """

        print('PROVENANCE FOR: {}'.format(fname))
        G = copy.deepcopy(self)
        
        f_node = G.FileNode(os.path.realpath(fname), Hasher()(fname), os.stat(fname)[stat.ST_MTIME])
        print('f_node=', f_node)

        hash2files, realpath2files = G._compute_hash_and_realpath_indices()
        print('hash2files:', '\n'.join(map(str, hash2files.get(f_node.file_hash, []))))
        print('realpath2files:', '\n'.join(map(str, realpath2files.get(f_node.realpath, []))))

        print('f_node in G? ', f_node in G)
        if f_node in G:
            print('in_degree is', G.in_degree[f_node])

        if f_node not in G:
            G.add_node(f_node, node_kind='file')

        if G.in_degree[f_node] < 1:
            print('trying to reconnect')
            f_node = G._reconstruct_fnode_maker(f_node)

        assert f_node in G

        assert G.in_degree[f_node] == 1, 'No provenance info for file {}'.format(fname)
        ancs = networkx.algorithms.dag.ancestors(G, f_node)
        for n in ancs:
            if G.is_snode(n):
                print(G.nodes[n]['step_record_fname'])

        if svgfile: G.write_svg(svgfile, nodes=list(ancs)+[f_node])
        
# end: class ProvenanceGraph(object)
