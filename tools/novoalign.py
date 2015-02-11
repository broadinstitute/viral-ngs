'''
    Novoalign aligner by Novocraft
    
    This is commercial software that has different licenses depending
    on use cases. As such, we do not have an auto-downloader. The user
    must have Novoalign pre-installed on their own and available
    either in $PATH or $NOVOALIGN_PATH.
'''

import tools, tools.picard, tools.samtools, util.file
import logging, os, os.path, subprocess, stat, gzip

log = logging.getLogger(__name__)

class NovoalignTool(tools.Tool) :
    def __init__(self, path=None):
        self.tool_version = None
        install_methods = []
        for novopath in [path, os.environ.get('NOVOALIGN_PATH'), '']:
            if novopath != None:
                install_methods.append(tools.PrexistingUnixCommand(
                    os.path.join(novopath, 'novoalign'),
                    require_executability=True))
        tools.Tool.__init__(self, install_methods = install_methods)
    
    def version(self):
        if self.tool_version==None:
            self._get_tool_version()
        return self.tool_version

    def _get_tool_version(self):
        tmpf = util.file.mkstempfname('.novohelp.txt')
        with open(tmpf, 'wt') as outf:
            subprocess.call([self.install_and_get_path()], stdout=outf)
        with open(tmpf, 'rt') as inf:
            self.tool_version = inf.readline().strip().split()[1]
        os.unlink(tmpf)
    
    def _fasta_to_idx_name(self, fasta):
        if not fasta.endswith('.fasta'):
            raise ValueError('input file %s must end with .fasta' % fasta)
        return fasta[:-6] + '.nix'
    
    def execute(self, inBam, refFasta, outBam,
        options=["-r", "Random"], min_qual=0, JVMmemory=None):
        ''' Execute Novoalign on BAM inputs and outputs.
            If the BAM contains multiple read groups, break up
            the input and perform Novoalign separately on each one
            (because Novoalign mangles read groups).
            Use Picard to sort and index the output BAM.
            If min_qual>0, use Samtools to filter on mapping quality.
        '''
        samtools = tools.samtools.SamtoolsTool()
        
        # Novoalign
        tmp_sam = util.file.mkstempfname('.novoalign.sam')
        cmd = [self.install_and_get_path(), '-f', inBam] + list(map(str, options))
        cmd = cmd + ['-F', 'BAMPE', '-d', self._fasta_to_idx_name(refFasta), '-o', 'SAM']
        log.debug(' '.join(cmd))
        with open(tmp_sam, 'wt') as outf:
            subprocess.check_call(cmd, stdout=outf)
        
        # Samtools filter (optional)
        if min_qual:
            tmp_bam2 = util.file.mkstempfname('.filtered.bam')
            cmd = [samtools.install_and_get_path(), 'view', '-b', '-S', '-1', '-q', str(min_qual), tmp_sam]
            log.debug('%s > %s' % (' '.join(cmd), tmp_bam2))
            with open(tmp_bam2, 'wb') as outf:
                subprocess.check_call(cmd, stdout=outf)
            os.unlink(tmp_sam)
            tmp_sam = tmp_bam2
        
        # Picard SortSam
        sorter = tools.picard.SortSamTool()
        sorter.execute(tmp_sam, outBam, sort_order='coordinate',
            picardOptions=['CREATE_INDEX=true', 'VALIDATION_STRINGENCY=SILENT'],
            JVMmemory=JVMmemory)
        

    def index_fasta(self, refFasta):
        ''' Index a FASTA file (reference genome) for use with Novoalign.
            The input file name must end in ".fasta". This will create a
            new ".nix" file in the same directory. If it already exists,
            it will be deleted and regenerated.
        '''
        novoindex = os.path.join(os.path.dirname(self.install_and_get_path()), 'novoindex')
        outfname = self._fasta_to_idx_name(refFasta)
        if os.path.isfile(outfname):
            os.unlink(outfname)
        cmd = [novoindex, outfname, refFasta]
        log.debug(' '.join(cmd))
        subprocess.check_call(cmd)
        try:
            mode = os.stat(outfname).st_mode & ~stat.S_IXUSR & ~stat.S_IXGRP & ~stat.S_IXOTH
            os.chmod(outfname, mode)
        except PermissionError:
            pass


class InvalidBamHeaderError(ValueError):
    pass
