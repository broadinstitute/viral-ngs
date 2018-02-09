"""Analysis of recorded metadata"""

import collections
import os
import os.path
import json
import fnmatch

import util.file
import util.misc
import functools
import itertools
import operator

from . import metadata_dir
from .file_arg import FileArg
from .hashing import Hasher
from .provenance_graph import ProvenanceGraph

import networkx
import networkx.algorithms.dag

##########################
# * Analysis of metadata 
##########################


# ** class Comp
class Comp(object):
    """"A Comp is a particular computation, represented by a subgraph of the ProvenanceGraph.
    For example, the steps needed to assemble a viral genome from a particular sample, and to compute metrics for the assembly.

    Fields:
       nodes: ProvenanceGraph nodes (both fnodes and snodes) comprising this computation
       main_inputs: fnode(s) denoting the main inputs to the computation, the data we're trying to analyze.
         For an assembly, the main inputs might be the raw reads files.  By contrast, things like depletion databases
         are more like parameters of the computation, rather than its inputs.
       main_outputs: fnode(s) denoting the main outputs of the computation.  This is in contrast to files that store metrics,
         log files, etc.
    """

    def __init__(self, nodes, main_inputs, main_outputs):
        self.nodes, self.main_inputs, self.main_outputs = nodes, main_inputs, main_outputs

    def __str__(self):
        return('Comp(main_inputs={}, main_outputs={})'.format(self.main_inputs, self.main_outputs))

# end: class Comp(object):

# class CompExtractor(object):
#     """Abstract base class for classes that extract particular kinds of Comps from the provenance graph."""
    
#     __metaclass__ = abc.ABCMeta

#     @abc.abstractmethod
#     def 


def compute_comp_attrs(G, comp):
    """Compute the attributes of a Comp, and gather them together into a single set of attribute-value pairs.  
    From each step, we gather the parameters of that step.

    Returns:
       the set of attribute-value pairs for the Comp `comp`
    """

    attrs={}
    for n in comp.nodes:
        if G.is_snode(n):
            na = G.nodes[n]
            step_name = na['step_name']  # might need to prepend cmd_module for uniqueness?

            def gather_files(val):
                """Return a flattened list of files denoted by this command argument (or empty list if it does not denote files)."""
                return (FileArg.is_from_dict(val) and [val]) \
                    or (isinstance(val, (list, tuple)) and functools.reduce(operator.concat, list(map(gather_files, val)), [])) or []

            for a, v in na['args'].items():
                if gather_files(v): continue
                if a in 'tmp_dir tmp_dirKeep loglevel'.split(): continue
                if step_name == 'impute_from_reference' and a == 'newName': continue
                attrs[step_name+'.'+a] = str(v)

            for a, v in na['metadata_from_cmd_return'].items():
                if a == '__metadata__': continue
                attrs[step_name+'.'+a] = str(v)

    return frozenset(attrs.items())

# ** extract_comps_assembly
def extract_comps_assembly(G):
    """From the provenance graph, extract Comps representing the assembly of one viral sample.

    Returns:
       a list of extracted comps, represented as Comp objects
    """

    extracted_comps = []

    beg_nodes = set([f for f in G.file_nodes if fnmatch.fnmatch(f.realpath, '*/data/00_raw/*.bam')])
    end_nodes = [f for f in G.file_nodes if fnmatch.fnmatch(f.realpath, '*/data/02_assembly/*.fasta')]
    for end_node in end_nodes:
        ancs = networkx.algorithms.dag.ancestors(G, end_node)
        ancs.add(end_node)
        ancs = set(ancs)
        beg = beg_nodes & ancs
        if beg:
            assert len(beg) == 1

            # add nodes that compute metrics
            final_assembly_metrics = [s for s in list(G.succ[end_node]) if G.nodes[s]['step_name']=='assembly_metrics']
            if final_assembly_metrics:
                final_assembly_metrics = sorted(final_assembly_metrics, key=lambda s: G.nodes[s]['run_info']['beg_time'], reverse=True)
                ancs.add(final_assembly_metrics[0])
                extracted_comps.append(Comp(nodes=ancs, main_inputs=list(beg), main_outputs=[end_node]))

    return extracted_comps

def group_comps_by_main_input(comps):
    """Group Comps by the contents of their main inputs.  This way we identify groups of comps where within each group,
    the comps all start with the same inputs."""

    def comp_main_inputs_contents(comp): return tuple(map(operator.attrgetter('file_hash'), comp.main_inputs))
    return [tuple(g) for k, g in itertools.groupby(sorted(comps, key=comp_main_inputs_contents), key=comp_main_inputs_contents)]

def report_comps_groups():

    G = ProvenanceGraph()
    G.load()

    grps = group_comps_by_main_input(extract_comps_assembly(G))
    print(sorted(map(len, grps)))
    for grp_num, g in enumerate(grps):
        g = list(g)
        if len(g) > 1:
            print('==============grp:')
            print('\n'.join(map(str, g)))
            for comp_num, comp in enumerate(g):
                G.write_svg('prov{}_{}.svg'.format(grp_num, comp_num), nodes=comp.nodes)
                if comp_num > 0:
                    print('symdiff:\n', '\n'.join(map(str, 
                                                    sorted(compute_comp_attrs(G, g[0]) ^ compute_comp_attrs(G, comp)))))


# end: def extract_paths():
