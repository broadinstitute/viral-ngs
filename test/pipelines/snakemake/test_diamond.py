#!/usr/bin/env python

import os
import sys

import pytest

import tools
from test.pipelines.snakemake import snake
from test.integration.test_diamond import * # for pytest fixtures

@pytest.mark.skipif(tools.is_osx(), reason="not currently tested under OSX")
@pytest.mark.skipif(sys.version_info < (3, 5), reason="Python version is too old for snakemake.")
def test_pipes(tmpdir_function, diamond_db, taxonomy_db, krona_db, input_bam):
    runner = snake.SnakemakeRunner(workdir=tmpdir_function)
    override_config = {
        'diamond_db': diamond_db,
        'taxonomy_db': taxonomy_db,
        'krona_db': krona_db,
    }
    runner.set_override_config(override_config)
    runner.setup()
    runner.link_samples([input_bam], destination='per_sample', link_transform=snake.rename_raw_bam)
    runner.create_sample_files(sample_files=['samples_metagenomics'])

    krona_out = join(runner.config['data_dir'], runner.config['subdirs']['metagenomics'],
                         '.'.join([os.path.splitext(os.path.basename(input_bam))[0], 'raw', 'diamond.krona.html']))

    diamond_out = join(
        runner.config['data_dir'], runner.config['subdirs']['metagenomics'],
        '.'.join([os.path.splitext(os.path.basename(input_bam))[0], 'raw', 'diamond.report'])
    )
    runner.run([krona_out])
    assert os.path.getsize(os.path.join(runner.workdir, diamond_out)) > 0
    assert os.path.getsize(os.path.join(runner.workdir, krona_out)) > 0