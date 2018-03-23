import assembly
import util.file
import Bio.SeqIO

@util.cmd_plugins.cmd_hookimpl
def cmd_compute_metadata_from_file_contents(fname):
    """Gather stats for a .fasta file"""

    print('FASTA PLUGIN: ', fname)

    result = {}
    if not util.file.ispipe(fname) and any(fname.endswith(ext) for ext in '.fasta .fa .fasta.gz .fa.gz'.split()):

        # should avoid reading whole into memory

        with util.file.open_or_gzopen(fname) as f:
            seqs = list(Bio.SeqIO.parse(f, 'fasta'))
        len_unambig = sum(map(assembly.unambig_count, seqs))
        len_tot = sum(map(len, seqs))
        result.update(fasta_n_seqs=len(seqs), fasta_len_unambig=len_unambig, fasta_len_tot=len_tot)

    return result
