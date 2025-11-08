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

    def execute(    self,
                    n_fastqs,
                    threads,
                    config_file,
                    keep_file,
                    unassigned_r1,
                    unassigned_r2,
                    summary_stats,
                    r1,
                    r2,
                    splitcode_opts=None,
                    gzip_output=False,
                    keep_r1_r2_suffixes=True,
                    predemux_r1_trim_5prime_num_bp=None,
                    predemux_r1_trim_3prime_num_bp=None,
                    predemux_r2_trim_5prime_num_bp=None,
                    predemux_r2_trim_3prime_num_bp=None
                ):
        """
        Execute the splitcode command with the provided parameters.

        :param n_fastqs: Number of FASTQ files (e.g., 2 for paired-end reads).
        :param threads: Number of threads to use.
        :param config_file: Path to the configuration file.
        :param keep_file: Path to the keep file.
        :param unassigned_r1: Path for unassigned R1 reads.
        :param unassigned_r2: Path for unassigned R2 reads.
        :param summary_stats: Path to the summary statistics file.
        :param r1: Input FASTQ file for R1 reads.
        :param r2: Input FASTQ file for R2 reads.
        :param splitcode_opts: additional parameters to pass to splitcode
        """
        splitcode_opts = splitcode_opts or []

        # documentation: 
        #   https://splitcode.readthedocs.io/en/latest/reference_guide.html#command-line-config-optional

        tool_cmd = [
            self.install_and_get_path(),
            f"--nFastqs={n_fastqs}",
            f"-t", str(threads),
            f"-c", config_file,
            f"--keep={keep_file}",
            f"--unassigned={unassigned_r1},{unassigned_r2}",
            f"--summary", summary_stats
        ]
        if gzip_output and "--gzip" not in splitcode_opts:
            tool_cmd.append("--gzip")
        if keep_r1_r2_suffixes and "--keep-r1-r2" not in splitcode_opts:
            tool_cmd.append("--keep-r1-r2")
        if predemux_r1_trim_3prime_num_bp or predemux_r2_trim_3prime_num_bp:
            r1_3prime_bp_to_remove = predemux_r1_trim_3prime_num_bp or 0
            r2_3prime_bp_to_remove = predemux_r2_trim_3prime_num_bp or 0
            tool_cmd.extend(["--trim-3",
                                f"{r1_3prime_bp_to_remove},{r2_3prime_bp_to_remove}"])
        if predemux_r1_trim_5prime_num_bp or predemux_r2_trim_5prime_num_bp:
            r1_5prime_bp_to_remove = predemux_r1_trim_5prime_num_bp or 0
            r2_5prime_bp_to_remove = predemux_r2_trim_5prime_num_bp or 0
            tool_cmd.extend(["--trim-5", 
                                f"{r1_5prime_bp_to_remove},{r2_5prime_bp_to_remove}"])

        tool_cmd.extend(splitcode_opts)
        tool_cmd.extend([r1,r2])

        log.debug('Running splitcode with command: %s', ' '.join(tool_cmd))
        try:
            return subprocess.check_call(tool_cmd)
        except subprocess.CalledProcessError as e:
            log.error(f'Splitcode failed with return code {e.returncode}')
            log.error(f'Command was: {" ".join(tool_cmd)}')
            raise

    def check_installation(self):
        """Ensure the tool is properly installed."""
        if not shutil.which(TOOL_NAME):
            raise FileNotFoundError(f"{TOOL_NAME} is not installed or not in PATH.")

    def run_splitcode(self, **kwargs):
        """
        Wrapper method to execute the splitcode command with named parameters.

        Expected kwargs:
        - n_fastqs, threads, config_file, keep_file, unassigned_r1, unassigned_r2, summary_stats, r1, r2
        """
        self.execute(
            kwargs.get('n_fastqs', 2),
            kwargs.get('threads', 1),
            kwargs['config_file'],
            kwargs['keep_file'],
            kwargs['unassigned_r1'],
            kwargs['unassigned_r2'],
            kwargs['summary_stats'],
            kwargs['r1'],
            kwargs['r2']
        )
