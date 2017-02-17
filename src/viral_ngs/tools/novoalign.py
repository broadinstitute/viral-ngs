'''
    Novoalign aligner by Novocraft

    This is commercial software that has different licenses depending
    on use cases. As such, we do not have an auto-downloader. The user
    must have Novoalign pre-installed on their own and available
    either in $PATH or $NOVOALIGN_PATH.
'''

import tools
import tools.picard
import tools.samtools
import util.file
import util.misc

import logging
import os
import os.path
import subprocess
import stat
import sys

_log = logging.getLogger(__name__)

TOOL_NAME = "novoalign"
TOOL_VERSION = "3.06.05"


class NovoalignTool(tools.Tool):

    def __init__(self, path=None, license_path=None):
        self.tool_version = None
        install_methods = []
        for novopath in [path, os.environ.get('NOVOALIGN_PATH'), '']:
            if novopath is not None:
                install_methods.append(
                    tools.PrexistingUnixCommand(
                        os.path.join(novopath, 'novoalign'),
                        require_executability=True
                    )
                )

        post_verify_command = None
        for novo_license_path in [license_path, os.environ.get("NOVOALIGN_LICENSE_PATH"), '']:
            if novo_license_path is not None and os.path.isfile(novo_license_path):
                # called relative to the conda bin/ directory
                uname = os.uname()
                # we ideally want the "update" copy operation
                # but only GNU cp has it. On OSX, the license will be copied each time.
                if uname[0] == 'Darwin':
                    copy_operation = ''
                else:
                    copy_operation = '-u'

                post_verify_command = "cp {copy_operation} {lic_path} ./".format(copy_operation=copy_operation, lic_path=novo_license_path)
                break # if we've found a license file, stop looking
            
        install_methods.append(tools.CondaPackage(TOOL_NAME, version=TOOL_VERSION, post_verify_command=post_verify_command))
        tools.Tool.__init__(self, install_methods=install_methods)

    def version(self):
        if self.tool_version is None:
            self._get_tool_version()
        return self.tool_version

    def _get_tool_version(self):
        tmpf = util.file.mkstempfname('.novohelp.txt')
        with open(tmpf, 'wt') as outf:
            util.misc.run_and_save([self.install_and_get_path(), "-V"], outf=outf, check=False)
        with open(tmpf, 'rt') as inf:
            self.tool_version = inf.readline().strip().split()[1]
        os.unlink(tmpf)

    def _fasta_to_idx_name(self, fasta):
        if not fasta.endswith('.fasta'):
            raise ValueError('input file %s must end with .fasta' % fasta)
        return fasta[:-6] + '.nix'

    def execute(self, inBam, refFasta, outBam, options=None, min_qual=0, JVMmemory=None):    # pylint: disable=W0221
        ''' Execute Novoalign on BAM inputs and outputs.
            If the BAM contains multiple read groups, break up
            the input and perform Novoalign separately on each one
            (because Novoalign mangles read groups).
            Use Picard to sort and index the output BAM.
            If min_qual>0, use Samtools to filter on mapping quality.
        '''
        options = options or ["-r", "Random"]

        samtools = tools.samtools.SamtoolsTool()

        # fetch list of RGs
        rgs = list(samtools.getReadGroups(inBam).keys())

        if len(rgs) == 0:
            # Can't do this
            raise InvalidBamHeaderError("{} lacks read groups".format(inBam))

        elif len(rgs) == 1:
            # Only one RG, keep it simple
            self.align_one_rg_bam(inBam, refFasta, outBam, options=options, min_qual=min_qual, JVMmemory=JVMmemory)

        else:
            # Multiple RGs, align one at a time and merge
            align_bams = []
            for rg in rgs:
                tmp_bam = util.file.mkstempfname('.{}.bam'.format(rg))
                self.align_one_rg_bam(
                    inBam,
                    refFasta,
                    tmp_bam,
                    rgid=rg,
                    options=options,
                    min_qual=min_qual,
                    JVMmemory=JVMmemory
                )
                if os.path.getsize(tmp_bam) > 0:
                    align_bams.append(tmp_bam)

            # Merge BAMs, sort, and index
            tools.picard.MergeSamFilesTool().execute(
                align_bams,
                outBam,
                picardOptions=['SORT_ORDER=coordinate', 'USE_THREADING=true', 'CREATE_INDEX=true'],
                JVMmemory=JVMmemory
            )
            for bam in align_bams:
                os.unlink(bam)

    def align_one_rg_bam(self, inBam, refFasta, outBam, rgid=None, options=None, min_qual=0, JVMmemory=None):
        ''' Execute Novoalign on BAM inputs and outputs.
            Requires that only one RG exists (will error otherwise).
            Use Picard to sort and index the output BAM.
            If min_qual>0, use Samtools to filter on mapping quality.
        '''
        options = options or ["-r", "Random"]

        samtools = tools.samtools.SamtoolsTool()

        # Require exactly one RG
        rgs = samtools.getReadGroups(inBam)
        if len(rgs) == 0:
            raise InvalidBamHeaderError("{} lacks read groups".format(inBam))
        elif len(rgs) == 1:
            if not rgid:
                rgid = list(rgs.keys())[0]
        elif not rgid:
            raise InvalidBamHeaderError("{} has {} read groups, but we require exactly one".format(inBam, len(rgs)))
        if rgid not in rgs:
            raise InvalidBamHeaderError("{} has read groups, but not {}".format(inBam, rgid))
        #rg = rgs[rgid]

        # Strip inBam to just one RG (if necessary)
        if len(rgs) == 1:
            one_rg_inBam = inBam
        else:
            # strip inBam to one read group
            tmp_bam = util.file.mkstempfname('.onebam.bam')
            samtools.view(['-b', '-r', rgid], inBam, tmp_bam)
            # special exit if this file is empty
            if samtools.count(tmp_bam) == 0:
                return
            # simplify BAM header otherwise Novoalign gets confused
            one_rg_inBam = util.file.mkstempfname('.{}.in.bam'.format(rgid))
            headerFile = util.file.mkstempfname('.{}.header.txt'.format(rgid))
            with open(headerFile, 'wt') as outf:
                for row in samtools.getHeader(inBam):
                    if len(row) > 0 and row[0] == '@RG':
                        if rgid != list(x[3:] for x in row if x.startswith('ID:'))[0]:
                            # skip all read groups that are not rgid
                            continue
                    outf.write('\t'.join(row) + '\n')
            samtools.reheader(tmp_bam, headerFile, one_rg_inBam)
            os.unlink(tmp_bam)
            os.unlink(headerFile)

        # Novoalign
        tmp_sam = util.file.mkstempfname('.novoalign.sam')
        tmp_sam_err = util.file.mkstempfname('.novoalign.sam.err')
        cmd = [self.install_and_get_path(), '-f', one_rg_inBam] + list(map(str, options))
        cmd = cmd + ['-F', 'BAM', '-d', self._fasta_to_idx_name(refFasta), '-o', 'SAM']
        _log.debug(' '.join(cmd))
        with open(tmp_sam, 'wt') as outf:
            util.misc.run_and_save(cmd, outf=outf)

        # Samtools filter (optional)
        if min_qual:
            tmp_bam2 = util.file.mkstempfname('.filtered.bam')
            cmd = [samtools.install_and_get_path(), 'view', '-b', '-S', '-1', '-q', str(min_qual), tmp_sam]
            _log.debug('%s > %s', ' '.join(cmd), tmp_bam2)
            with open(tmp_bam2, 'wb') as outf:
                util.misc.run_and_save(cmd, outf=outf)
            os.unlink(tmp_sam)
            tmp_sam = tmp_bam2

        # Picard SortSam
        sorter = tools.picard.SortSamTool()
        sorter.execute(
            tmp_sam,
            outBam,
            sort_order='coordinate',
            picardOptions=['CREATE_INDEX=true', 'VALIDATION_STRINGENCY=SILENT'],
            JVMmemory=JVMmemory
        )

    def index_fasta(self, refFasta, k=None, s=None):
        ''' Index a FASTA file (reference genome) for use with Novoalign.
            The input file name must end in ".fasta". This will create a
            new ".nix" file in the same directory. If it already exists,
            it will be deleted and regenerated.
        '''
        novoindex = os.path.join(os.path.dirname(self.install_and_get_path()), 'novoindex')
        outfname = self._fasta_to_idx_name(refFasta)
        if os.path.isfile(outfname):
            os.unlink(outfname)
        cmd = [novoindex]
        if k is not None:
            cmd.extend(['-k', str(k)])
        if s is not None:
            cmd.extend(['-s', str(s)])
        cmd.extend([outfname, refFasta])
        _log.debug(' '.join(cmd))
        subprocess.check_call(cmd)
        try:
            mode = os.stat(outfname).st_mode & ~stat.S_IXUSR & ~stat.S_IXGRP & ~stat.S_IXOTH
            os.chmod(outfname, mode)
        except (IOError, OSError):
            _log.warning('could not chmod "%s", this is likely OK', refFasta)


class InvalidBamHeaderError(ValueError):
    pass
