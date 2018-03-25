import tools.samtools
import util.file
import util.cmd_plugins

@util.cmd_plugins.cmd_hookimpl
def cmd_compute_metadata_from_file_contents(fname):
    """Gather stats for a .bam file"""

    result = {}
    if fname.endswith('.bam') and not util.file.ispipe(fname):
        result['bam_count'] = tools.samtools.SamtoolsTool().count(fname)
    return result
