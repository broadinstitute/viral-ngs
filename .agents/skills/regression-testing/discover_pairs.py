#!/usr/bin/env python3
"""Discover comparable old/new sample pairs by crawling GCS Cromwell output directories.

For each submission, finds assembly_metadata TSV files named
``assembly_metadata-<sample>.tsv``, extracts sample names from filenames,
and outputs the intersection as a JSON mapping.

Usage:
    python discover_pairs.py \
        --bucket fc-XXXXXXXX-... \
        --old-sub <old-submission-id> \
        --new-sub <new-submission-id> \
        -o pairs.json
"""
import argparse
import json
import logging
import re
import subprocess
import sys
from collections import defaultdict

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')
log = logging.getLogger(__name__)


def gcloud_ls(path):
    """List GCS path, return list of URIs. Returns [] on error."""
    try:
        result = subprocess.run(
            ['gcloud', 'storage', 'ls', path],
            capture_output=True, text=True, timeout=60
        )
        if result.returncode != 0:
            return []
        return [line.strip() for line in result.stdout.strip().split('\n') if line.strip()]
    except Exception as e:
        log.warning(f"gcloud ls failed for {path}: {e}")
        return []


def find_tsv_in_call_dir(call_dir_uri):
    """Find assembly_metadata TSV in a call directory, handling attempt-N subdirs.

    Returns (sample_name, tsv_uri) or (None, None).
    """
    items = gcloud_ls(call_dir_uri)

    # Check for attempt-N subdirectories
    def attempt_sort_key(path):
        match = re.search(r'/attempt-(\d+)', path)
        return int(match.group(1)) if match else 0
    attempt_dirs = sorted([i for i in items if '/attempt-' in i],
                          key=attempt_sort_key, reverse=True)
    tsv_files = [i for i in items if i.endswith('.tsv')]

    # If there are attempt dirs, check the highest attempt first
    if attempt_dirs:
        for attempt_dir in attempt_dirs:
            attempt_items = gcloud_ls(attempt_dir)
            attempt_tsvs = [i for i in attempt_items if i.endswith('.tsv')]
            if attempt_tsvs:
                tsv_files = attempt_tsvs
                break

    for tsv in tsv_files:
        match = re.search(r'assembly_metadata-(.+)\.tsv$', tsv)
        if match:
            return match.group(1), tsv

    return None, None


def discover_submission_tsvs(bucket, submission_id):
    """Find all assembly_stats TSVs for a submission.

    Returns dict: sample_name -> tsv_gcs_uri
    """
    base = f"gs://{bucket}/submissions/{submission_id}/assemble_denovo_metagenomic/"
    log.info(f"Listing workflow directories in {base}")

    wf_dirs = gcloud_ls(base)
    log.info(f"Found {len(wf_dirs)} workflow directories")

    results = {}
    for i, wf_dir in enumerate(wf_dirs):
        if i % 20 == 0:
            log.info(f"  Scanning workflow {i+1}/{len(wf_dirs)}...")

        for call_name in ['call-assembly_stats_non_empty', 'call-assembly_stats_empty']:
            call_dir = f"{wf_dir}{call_name}/"
            sample_name, tsv_uri = find_tsv_in_call_dir(call_dir)
            if sample_name:
                if sample_name in results:
                    log.warning(f"Duplicate sample {sample_name} -- keeping first occurrence")
                else:
                    results[sample_name] = tsv_uri
                break

    log.info(f"Found TSVs for {len(results)} samples in submission {submission_id[:8]}")
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--bucket', required=True,
                        help='Terra workspace GCS bucket ID (e.g., fc-XXXXXXXX-...)')
    parser.add_argument('--old-sub', required=True,
                        help='Old submission ID (main branch)')
    parser.add_argument('--new-sub', required=True,
                        help='New submission ID (feature branch)')
    parser.add_argument('--output', '-o', required=True,
                        help='Output JSON file path')
    args = parser.parse_args()

    log.info(f"Old submission: {args.old_sub[:8]}")
    log.info(f"New submission: {args.new_sub[:8]}")

    old_tsvs = discover_submission_tsvs(args.bucket, args.old_sub)
    new_tsvs = discover_submission_tsvs(args.bucket, args.new_sub)

    # Find intersection
    common_samples = sorted(set(old_tsvs.keys()) & set(new_tsvs.keys()))
    old_only = sorted(set(old_tsvs.keys()) - set(new_tsvs.keys()))
    new_only = sorted(set(new_tsvs.keys()) - set(old_tsvs.keys()))

    log.info(f"Old-only samples: {len(old_only)}")
    log.info(f"New-only samples: {len(new_only)}")
    log.info(f"Intersecting samples: {len(common_samples)}")

    if old_only:
        log.info(f"  Old-only: {old_only[:5]}{'...' if len(old_only) > 5 else ''}")
    if new_only:
        log.info(f"  New-only: {new_only[:5]}{'...' if len(new_only) > 5 else ''}")

    pairs = {}
    for sample in common_samples:
        pairs[sample] = {
            'old_tsv': old_tsvs[sample],
            'new_tsv': new_tsvs[sample],
        }

    output = {
        'bucket': args.bucket,
        'old_submission': args.old_sub,
        'new_submission': args.new_sub,
        'old_sample_count': len(old_tsvs),
        'new_sample_count': len(new_tsvs),
        'paired_count': len(pairs),
        'old_only': old_only,
        'new_only': new_only,
        'pairs': pairs,
    }

    with open(args.output, 'w') as f:
        json.dump(output, f, indent=2)

    log.info(f"Wrote {len(pairs)} pairs to {args.output}")


if __name__ == '__main__':
    main()
