import logging
import shutil
import subprocess
import tools

TOOL_NAME = 'splitcode'

log = logging.getLogger(__name__)

class SplitCodeTool(tools.Tool):

    def __init__(self, install_methods=None):
        if install_methods is None:
            install_methods = [tools.PrexistingUnixCommand(shutil.which(TOOL_NAME), require_executability=True)]
        super(SplitCodeTool, self).__init__(install_methods=install_methods)

    def _get_tool_version(self):
        self.tool_version = subprocess.check_output([self.install_and_get_path(), '--version']).decode('UTF-8').strip()

    def execute(self, nFastqs, threads, config_file, keep_file, unassigned_r1, unassigned_r2, summary_stats, r1, r2):
        """
        Execute the splitcode command with the provided parameters.

        :param nFastqs: Number of FASTQ files (e.g., 2 for paired-end reads).
        :param threads: Number of threads to use.
        :param config_file: Path to the configuration file.
        :param keep_file: Path to the keep file.
        :param unassigned_r1: Path for unassigned R1 reads.
        :param unassigned_r2: Path for unassigned R2 reads.
        :param summary_stats: Path to the summary statistics file.
        :param r1: Input FASTQ file for R1 reads.
        :param r2: Input FASTQ file for R2 reads.
        """
        tool_cmd = [
            self.install_and_get_path(),
            f"--nFastqs={nFastqs}",
            f"-t", str(threads),
            f"-c", config_file,
            f"--keep={keep_file}",
            f"--unassigned={unassigned_r1},{unassigned_r2}",
            "--no-output",
            "--no-outb",
            f"--summary", summary_stats,
            r1,
            r2,
        ]

        log.debug('Running splitcode with command: %s', ' '.join(tool_cmd))
        subprocess.check_call(tool_cmd)

    def check_installation(self):
        """Ensure the tool is properly installed."""
        if not shutil.which(TOOL_NAME):
            raise FileNotFoundError(f"{TOOL_NAME} is not installed or not in PATH.")

    def run_splitcode(self, **kwargs):
        """
        Wrapper method to execute the splitcode command with named parameters.

        Expected kwargs:
        - nFastqs, threads, config_file, keep_file, unassigned_r1, unassigned_r2, summary_stats, r1, r2
        """
        self.execute(
            kwargs.get('nFastqs', 2),
            kwargs.get('threads', 1),
            kwargs['config_file'],
            kwargs['keep_file'],
            kwargs['unassigned_r1'],
            kwargs['unassigned_r2'],
            kwargs['summary_stats'],
            kwargs['r1'],
            kwargs['r2']
        )
