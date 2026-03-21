#!/usr/bin/env python3
"""Compare assembly outputs between old and new code for a single sample pair.

Takes two GCS URIs pointing at assembly_metadata TSV files (old and new),
downloads them, compares metrics, aligns FASTAs with mafft, and outputs a JSON result.
"""
import argparse
import csv
import io
import json
import logging
import os
import subprocess

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')
log = logging.getLogger(__name__)

METRICS_COLS = [
    'assembly_length', 'assembly_length_unambiguous', 'reads_aligned',
    'mean_coverage', 'percent_reference_covered', 'reference_length',
    'scaffolding_num_segments_recovered', 'reference_num_segments_required',
]
FLOAT_COLS = {'mean_coverage', 'percent_reference_covered'}
INT_COLS = {'assembly_length', 'assembly_length_unambiguous', 'reads_aligned',
            'reference_length', 'scaffolding_num_segments_recovered',
            'reference_num_segments_required'}


def gcloud_cat(gcs_uri):
    """Read a GCS file's contents as a string."""
    result = subprocess.run(
        ['gcloud', 'storage', 'cat', gcs_uri],
        capture_output=True, text=True, timeout=120
    )
    if result.returncode != 0:
        raise RuntimeError(f"gcloud storage cat failed for {gcs_uri}: {result.stderr.strip()}")
    return result.stdout


def gcloud_cp(gcs_uri, local_path):
    """Download a GCS file to local path."""
    result = subprocess.run(
        ['gcloud', 'storage', 'cp', gcs_uri, local_path],
        capture_output=True, text=True, timeout=120
    )
    if result.returncode != 0:
        raise RuntimeError(f"gcloud storage cp failed for {gcs_uri}: {result.stderr.strip()}")


def parse_tsv(content):
    """Parse assembly_metadata TSV content into a dict keyed by assembly_id.

    Returns dict: assembly_id -> {col: value, ...}
    """
    reader = csv.DictReader(io.StringIO(content), delimiter='\t')
    rows = {}
    for row in reader:
        # The first column is 'entity:assembly_id'
        assembly_id = row.get('entity:assembly_id', '').strip()
        if not assembly_id:
            continue
        # Parse numeric columns
        parsed = dict(row)
        for col in FLOAT_COLS:
            if col in parsed and parsed[col]:
                try:
                    parsed[col] = float(parsed[col])
                except (ValueError, TypeError):
                    parsed[col] = None
        for col in INT_COLS:
            if col in parsed and parsed[col]:
                try:
                    parsed[col] = int(parsed[col])
                except (ValueError, TypeError):
                    parsed[col] = None
        rows[assembly_id] = parsed
    return rows


def parse_fasta(content):
    """Parse a FASTA string into list of (header, sequence) tuples."""
    sequences = []
    current_header = None
    current_seq = []
    for line in content.strip().split('\n'):
        line = line.strip()
        if not line:
            continue
        if line.startswith('>'):
            if current_header is not None:
                sequences.append((current_header, ''.join(current_seq)))
            current_header = line[1:].strip()
            current_seq = []
        else:
            current_seq.append(line)
    if current_header is not None:
        sequences.append((current_header, ''.join(current_seq)))
    return sequences


def run_mafft_pair(old_seq, new_seq, work_dir, pair_id='0'):
    """Run mafft on a single pair of sequences. Returns aligned (old_seq, new_seq) strings."""
    combined = os.path.join(work_dir, f'combined_{pair_id}.fasta')
    with open(combined, 'w') as f:
        f.write(f'>old_{pair_id}\n{old_seq}\n>new_{pair_id}\n{new_seq}\n')

    try:
        result = subprocess.run(
            ['mafft', '--auto', '--preservecase', '--quiet', '--thread', '1', combined],
            capture_output=True, text=True, timeout=300
        )
    finally:
        os.unlink(combined)
    if result.returncode != 0:
        raise RuntimeError(f"mafft failed: {result.stderr.strip()}")

    seqs = parse_fasta(result.stdout)
    if len(seqs) != 2:
        raise RuntimeError(f"Expected 2 sequences from mafft, got {len(seqs)}")
    return seqs[0][1], seqs[1][1]


def align_and_analyze_fastas(old_fasta_path, new_fasta_path, work_dir):
    """Align old vs new FASTAs, handling multi-segment genomes.

    For single-segment: runs one mafft alignment and analyzes it.
    For multi-segment: pairs segments by header, aligns each pair separately,
    analyzes each independently (so terminal effects stay terminal), and
    aggregates the stats.

    Returns (aligned_fasta_path, stats_dict).
    """
    with open(old_fasta_path) as f:
        old_fasta_text = f.read()
    with open(new_fasta_path) as f:
        new_fasta_text = f.read()
    old_seqs = parse_fasta(old_fasta_text)
    new_seqs = parse_fasta(new_fasta_text)

    if len(old_seqs) == 1 and len(new_seqs) == 1:
        # Simple case: single segment
        aln_old, aln_new = run_mafft_pair(old_seqs[0][1], new_seqs[0][1], work_dir)
        aligned = os.path.join(work_dir, 'aligned.fasta')
        with open(aligned, 'w') as f:
            f.write(f'>old_{old_seqs[0][0]}\n{aln_old}\n>new_{new_seqs[0][0]}\n{aln_new}\n')
        stats = analyze_alignment_seqs(aln_old, aln_new)
        return aligned, stats

    # Multi-segment: pair by header name
    old_dict = {h: s for h, s in old_seqs}
    new_dict = {h: s for h, s in new_seqs}
    common_headers = [h for h in old_dict if h in new_dict]

    if not common_headers:
        if len(old_seqs) == len(new_seqs):
            log.info(f"  Multi-segment: headers don't match, pairing by position ({len(old_seqs)} segments)")
            common_headers = None
        else:
            raise RuntimeError(
                f"Cannot pair segments: {len(old_seqs)} old vs {len(new_seqs)} new, no matching headers")

    # Align each pair and analyze independently
    segment_pairs = []
    if common_headers is not None:
        for i, h in enumerate(common_headers):
            aln_old, aln_new = run_mafft_pair(old_dict[h], new_dict[h], work_dir, pair_id=str(i))
            segment_pairs.append((h, aln_old, aln_new))
    else:
        for i in range(len(old_seqs)):
            aln_old, aln_new = run_mafft_pair(old_seqs[i][1], new_seqs[i][1], work_dir, pair_id=str(i))
            segment_pairs.append((old_seqs[i][0], aln_old, aln_new))

    n_segments = len(segment_pairs)
    log.info(f"  Multi-segment alignment: {n_segments} segments aligned")

    # Analyze each segment independently, then aggregate
    agg = {
        'alignment_length': 0, 'internal_length': 0,
        'matches': 0, 'snps': 0,
        'internal_insertions': 0, 'internal_deletions': 0,
        'internal_insertion_events': 0, 'internal_deletion_events': 0,
        'ambiguity_diffs': 0,
        'terminal_old_left': 0, 'terminal_new_left': 0,
        'terminal_old_right': 0, 'terminal_new_right': 0,
        'terminal_extensions_old': 0, 'terminal_extensions_new': 0,
        'terminal_extension_events_old': 0, 'terminal_extension_events_new': 0,
    }
    per_segment = []
    for header, aln_old, aln_new in segment_pairs:
        seg_stats = analyze_alignment_seqs(aln_old, aln_new)
        per_segment.append({'segment': header, **seg_stats})
        for key in agg:
            agg[key] += seg_stats.get(key, 0)

    # Compute aggregate identity
    total_bases_compared = agg['matches'] + agg['snps'] + agg['ambiguity_diffs']
    agg['identity'] = agg['matches'] / total_bases_compared if total_bases_compared > 0 else 1.0
    agg['n_segments'] = n_segments
    agg['per_segment'] = per_segment

    # Write per-segment alignment file (for review)
    aligned = os.path.join(work_dir, 'aligned.fasta')
    with open(aligned, 'w') as f:
        for header, aln_old, aln_new in segment_pairs:
            f.write(f'>old_{header}\n{aln_old}\n>new_{header}\n{aln_new}\n')

    return aligned, agg


def analyze_alignment_seqs(old_seq_str, new_seq_str):
    """Analyze a pairwise alignment given two aligned sequence strings.

    Returns dict with alignment statistics.
    """
    old_seq = old_seq_str.upper()
    new_seq = new_seq_str.upper()
    aln_len = len(old_seq)

    if len(new_seq) != aln_len:
        raise RuntimeError(f"Alignment sequences differ in length: {aln_len} vs {len(new_seq)}")

    ACGT = set('ACGT')

    # Find the internal region (where both sequences have bases)
    left_bound = 0
    while left_bound < aln_len and (old_seq[left_bound] == '-' or new_seq[left_bound] == '-'):
        left_bound += 1

    right_bound = aln_len - 1
    while right_bound >= 0 and (old_seq[right_bound] == '-' or new_seq[right_bound] == '-'):
        right_bound -= 1

    # Terminal extension counts (bp and events)
    terminal_old_left = 0
    terminal_new_left = 0
    terminal_old_right = 0
    terminal_new_right = 0
    # Track events: a contiguous run of gaps on one side = 1 event
    terminal_old_left_events = 0
    terminal_new_left_events = 0
    terminal_old_right_events = 0
    terminal_new_right_events = 0

    prev_old_gap = False
    prev_new_gap = False
    for i in range(left_bound):
        if old_seq[i] != '-' and new_seq[i] == '-':
            terminal_old_left += 1
            if not prev_old_gap:
                terminal_old_left_events += 1
            prev_old_gap = True
            prev_new_gap = False
        elif new_seq[i] != '-' and old_seq[i] == '-':
            terminal_new_left += 1
            if not prev_new_gap:
                terminal_new_left_events += 1
            prev_new_gap = True
            prev_old_gap = False
        else:
            prev_old_gap = False
            prev_new_gap = False

    prev_old_gap = False
    prev_new_gap = False
    for i in range(right_bound + 1, aln_len):
        if old_seq[i] != '-' and new_seq[i] == '-':
            terminal_old_right += 1
            if not prev_old_gap:
                terminal_old_right_events += 1
            prev_old_gap = True
            prev_new_gap = False
        elif new_seq[i] != '-' and old_seq[i] == '-':
            terminal_new_right += 1
            if not prev_new_gap:
                terminal_new_right_events += 1
            prev_new_gap = True
            prev_old_gap = False
        else:
            prev_old_gap = False
            prev_new_gap = False

    # Analyze internal region (bp counts and event counts)
    matches = 0
    snps = 0
    internal_insertions = 0  # bp count
    internal_deletions = 0   # bp count
    internal_insertion_events = 0
    internal_deletion_events = 0
    ambiguity_diffs = 0
    in_insertion = False
    in_deletion = False

    for i in range(left_bound, right_bound + 1):
        o = old_seq[i]
        n = new_seq[i]

        if o == '-' and n != '-':
            internal_insertions += 1
            if not in_insertion:
                internal_insertion_events += 1
                in_insertion = True
            in_deletion = False
        elif o != '-' and n == '-':
            internal_deletions += 1
            if not in_deletion:
                internal_deletion_events += 1
                in_deletion = True
            in_insertion = False
        else:
            in_insertion = False
            in_deletion = False
            if o == n:
                if o != '-':
                    matches += 1
            elif o in ACGT and n in ACGT:
                snps += 1
            else:
                ambiguity_diffs += 1

    total_internal = right_bound - left_bound + 1 if right_bound >= left_bound else 0
    total_bases_compared = matches + snps + ambiguity_diffs
    identity = matches / total_bases_compared if total_bases_compared > 0 else 1.0

    return {
        'alignment_length': aln_len,
        'internal_length': total_internal,
        'matches': matches,
        'snps': snps,
        'internal_insertions': internal_insertions,
        'internal_deletions': internal_deletions,
        'internal_insertion_events': internal_insertion_events,
        'internal_deletion_events': internal_deletion_events,
        'ambiguity_diffs': ambiguity_diffs,
        'terminal_old_left': terminal_old_left,
        'terminal_new_left': terminal_new_left,
        'terminal_old_right': terminal_old_right,
        'terminal_new_right': terminal_new_right,
        'terminal_extensions_old': terminal_old_left + terminal_old_right,
        'terminal_extensions_new': terminal_new_left + terminal_new_right,
        'terminal_extension_events_old': terminal_old_left_events + terminal_old_right_events,
        'terminal_extension_events_new': terminal_new_left_events + terminal_new_right_events,
        'identity': identity,
    }


def analyze_alignment(aligned_fasta_path):
    """Analyze pairwise alignment from mafft output file. Thin wrapper around analyze_alignment_seqs."""
    seqs = parse_fasta(open(aligned_fasta_path).read())
    if len(seqs) != 2:
        raise RuntimeError(f"Expected 2 sequences in alignment, got {len(seqs)}")
    return analyze_alignment_seqs(seqs[0][1], seqs[1][1])


def compare_assembly(assembly_id, old_row, new_row, work_dir):
    """Compare one assembly (old vs new). Downloads FASTAs, runs mafft, returns result dict."""
    result = {
        'assembly_id': assembly_id,
        'taxid': old_row.get('taxid', ''),
        'tax_name': old_row.get('tax_name', ''),
        'old_metrics': {},
        'new_metrics': {},
        'deltas': {},
        'alignment': None,
        'error': None,
    }

    # Extract metrics
    for col in METRICS_COLS:
        old_val = old_row.get(col)
        new_val = new_row.get(col)
        result['old_metrics'][col] = old_val
        result['new_metrics'][col] = new_val
        if old_val is not None and new_val is not None:
            try:
                result['deltas'][col] = float(new_val) - float(old_val)
            except (ValueError, TypeError):
                result['deltas'][col] = None

    # Get FASTA paths
    old_fasta_uri = old_row.get('assembly_fasta', '').strip()
    new_fasta_uri = new_row.get('assembly_fasta', '').strip()

    if not old_fasta_uri or not new_fasta_uri:
        log.info(f"  {assembly_id}: skipping alignment (missing FASTA URI)")
        return result

    # Download FASTAs
    aln_dir = os.path.join(work_dir, f'aln_{assembly_id}')
    os.makedirs(aln_dir, exist_ok=True)

    old_fasta_local = os.path.join(aln_dir, 'old.fasta')
    new_fasta_local = os.path.join(aln_dir, 'new.fasta')

    try:
        log.debug(f"  Downloading old FASTA: {old_fasta_uri}")
        gcloud_cp(old_fasta_uri, old_fasta_local)
        log.debug(f"  Downloading new FASTA: {new_fasta_uri}")
        gcloud_cp(new_fasta_uri, new_fasta_local)

        # Check if FASTAs are non-empty
        if os.path.getsize(old_fasta_local) == 0 or os.path.getsize(new_fasta_local) == 0:
            log.info(f"  {assembly_id}: skipping alignment (empty FASTA)")
            return result

        # Align and analyze (handles multi-segment genomes independently)
        log.debug(f"  Running alignment for {assembly_id}")
        aligned_path, alignment_stats = align_and_analyze_fastas(old_fasta_local, new_fasta_local, aln_dir)
        result['alignment'] = alignment_stats

        identity = alignment_stats['identity']
        snps = alignment_stats['snps']
        indels = alignment_stats['internal_insertions'] + alignment_stats['internal_deletions']
        log.info(f"  {assembly_id}: identity={identity:.6f}, snps={snps}, indels={indels}, "
                 f"terminal_ext_old={alignment_stats['terminal_extensions_old']}, "
                 f"terminal_ext_new={alignment_stats['terminal_extensions_new']}")

    except Exception as e:
        log.error(f"  {assembly_id}: alignment failed: {e}")
        result['error'] = str(e)
    finally:
        # Cleanup downloaded files
        for f in [old_fasta_local, new_fasta_local]:
            if os.path.exists(f):
                os.unlink(f)
        # Keep aligned file only if identity < 99.9%
        aligned_file = os.path.join(aln_dir, 'aligned.fasta')
        if os.path.exists(aligned_file):
            if result['alignment'] and result['alignment']['identity'] >= 0.999:
                os.unlink(aligned_file)
            else:
                log.info(f"  Keeping alignment file for review: {aligned_file}")
        # Remove empty dir
        try:
            os.rmdir(aln_dir)
        except OSError:
            pass  # dir not empty (kept alignment file)

    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--old-tsv', required=True, help='GCS URI of old assembly_metadata TSV')
    parser.add_argument('--new-tsv', required=True, help='GCS URI of new assembly_metadata TSV')
    parser.add_argument('--work-dir', required=True, help='Working directory for temp files')
    parser.add_argument('--output-json', required=True, help='Output JSON file path')
    parser.add_argument('--verbose', '-v', action='store_true', help='Enable debug logging')
    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    os.makedirs(args.work_dir, exist_ok=True)

    # Download and parse TSVs
    log.info(f"Downloading old TSV: {args.old_tsv}")
    old_content = gcloud_cat(args.old_tsv)
    old_rows = parse_tsv(old_content)

    log.info(f"Downloading new TSV: {args.new_tsv}")
    new_content = gcloud_cat(args.new_tsv)
    new_rows = parse_tsv(new_content)

    # Extract sample_id from first row (or from TSV filename)
    sample_id = 'unknown'
    if old_rows:
        first_row = next(iter(old_rows.values()))
        sample_id = first_row.get('sample_id', 'unknown')
    elif new_rows:
        first_row = next(iter(new_rows.values()))
        sample_id = first_row.get('sample_id', 'unknown')

    log.info(f"Sample: {sample_id}")
    log.info(f"Old assemblies: {len(old_rows)}, New assemblies: {len(new_rows)}")

    # Find intersecting assembly_ids
    common_ids = sorted(set(old_rows.keys()) & set(new_rows.keys()))
    old_only_ids = sorted(set(old_rows.keys()) - set(new_rows.keys()))
    new_only_ids = sorted(set(new_rows.keys()) - set(old_rows.keys()))

    if old_only_ids:
        log.info(f"  Assemblies only in old: {old_only_ids}")
    if new_only_ids:
        log.info(f"  Assemblies only in new: {new_only_ids}")
    log.info(f"  Assemblies in common: {len(common_ids)}")

    # Compare each intersecting assembly
    comparisons = []
    for aid in common_ids:
        comp = compare_assembly(aid, old_rows[aid], new_rows[aid], args.work_dir)
        comparisons.append(comp)

    # Build output
    output = {
        'sample_id': sample_id,
        'old_tsv_uri': args.old_tsv,
        'new_tsv_uri': args.new_tsv,
        'old_assembly_count': len(old_rows),
        'new_assembly_count': len(new_rows),
        'assembly_count_match': len(old_rows) == len(new_rows),
        'assemblies_only_in_old': old_only_ids,
        'assemblies_only_in_new': new_only_ids,
        'comparisons': comparisons,
    }

    with open(args.output_json, 'w') as f:
        json.dump(output, f, indent=2)

    log.info(f"Wrote results to {args.output_json}")


if __name__ == '__main__':
    main()
