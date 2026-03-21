#!/usr/bin/env python3
"""Generate regression testing report with plots from per-sample JSON results.

Aggregates all comparison results, produces summary TSV, plots, and markdown report.
"""
import argparse
import glob
import json
import logging
import os

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')
log = logging.getLogger(__name__)

# Delay imports so script can show usage without these deps
def get_deps():
    import pandas as pd
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    return pd, plt


def load_results(results_dir):
    """Load all per-sample JSON results."""
    results = []
    for path in sorted(glob.glob(os.path.join(results_dir, '*.json'))):
        try:
            with open(path) as f:
                results.append(json.load(f))
        except Exception as e:
            log.warning(f"Failed to load {path}: {e}")
    return results


def build_comparison_table(results, workspace_name):
    """Build a flat table of per-assembly comparisons."""
    rows = []
    for r in results:
        sample_id = r.get('sample_id', 'unknown')
        for comp in r.get('comparisons', []):
            row = {
                'workspace': workspace_name,
                'sample_id': sample_id,
                'assembly_id': comp.get('assembly_id', ''),
                'taxid': comp.get('taxid', ''),
                'tax_name': comp.get('tax_name', ''),
            }
            # Metrics
            for col in ['percent_reference_covered', 'mean_coverage',
                        'assembly_length_unambiguous', 'assembly_length',
                        'reads_aligned', 'reference_length']:
                old_val = comp.get('old_metrics', {}).get(col)
                new_val = comp.get('new_metrics', {}).get(col)
                delta = comp.get('deltas', {}).get(col)
                row[f'old_{col}'] = old_val
                row[f'new_{col}'] = new_val
                row[f'delta_{col}'] = delta

            # Alignment stats
            aln = comp.get('alignment')
            if aln:
                row['alignment_identity'] = aln.get('identity')
                row['snp_count'] = aln.get('snps', 0)
                row['internal_insertions'] = aln.get('internal_insertions', 0)
                row['internal_deletions'] = aln.get('internal_deletions', 0)
                row['indel_count_bp'] = aln.get('internal_insertions', 0) + aln.get('internal_deletions', 0)
                row['indel_count_events'] = aln.get('internal_insertion_events', 0) + aln.get('internal_deletion_events', 0)
                row['terminal_extensions_old'] = aln.get('terminal_extensions_old', 0)
                row['terminal_extensions_new'] = aln.get('terminal_extensions_new', 0)
                row['terminal_extension_events_old'] = aln.get('terminal_extension_events_old', 0)
                row['terminal_extension_events_new'] = aln.get('terminal_extension_events_new', 0)
                row['ambiguity_diffs'] = aln.get('ambiguity_diffs', 0)
            else:
                row['alignment_identity'] = None
                row['snp_count'] = None
                row['indel_count_bp'] = None
                row['indel_count_events'] = None
                row['terminal_extensions_old'] = None
                row['terminal_extensions_new'] = None
                row['terminal_extension_events_old'] = None
                row['terminal_extension_events_new'] = None
                row['ambiguity_diffs'] = None

            row['error'] = comp.get('error')
            rows.append(row)
    return rows


def build_sample_summary(results):
    """Build sample-level summary table."""
    rows = []
    for r in results:
        rows.append({
            'sample_id': r.get('sample_id', 'unknown'),
            'old_assembly_count': r.get('old_assembly_count', 0),
            'new_assembly_count': r.get('new_assembly_count', 0),
            'assembly_count_match': r.get('assembly_count_match', False),
            'assemblies_only_in_old': len(r.get('assemblies_only_in_old', [])),
            'assemblies_only_in_new': len(r.get('assemblies_only_in_new', [])),
            'num_comparisons': len(r.get('comparisons', [])),
        })
    return rows


def generate_plots(df, plot_dir):
    """Generate all plots from the comparison dataframe."""
    _, plt = get_deps()
    os.makedirs(plot_dir, exist_ok=True)

    # Filter to rows with actual assemblies
    df_asm = df[df['old_percent_reference_covered'].notna() &
                df['new_percent_reference_covered'].notna()].copy()

    if len(df_asm) == 0:
        log.warning("No assembly comparisons with metrics — skipping plots")
        return

    # 1. Percent reference covered scatter
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.scatter(df_asm['old_percent_reference_covered'] * 100,
               df_asm['new_percent_reference_covered'] * 100,
               alpha=0.5, s=20)
    lims = [0, 105]
    ax.plot(lims, lims, 'r--', alpha=0.5, label='y=x')
    ax.set_xlabel('Old % Reference Covered')
    ax.set_ylabel('New % Reference Covered')
    ax.set_title('Percent Reference Covered: Old vs New')
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(plot_dir, 'pct_ref_covered_scatter.png'), dpi=150)
    plt.close(fig)

    # 2. Mean coverage scatter (log scale)
    df_cov = df_asm[df_asm['old_mean_coverage'].notna() &
                     df_asm['new_mean_coverage'].notna() &
                     (df_asm['old_mean_coverage'] > 0) &
                     (df_asm['new_mean_coverage'] > 0)].copy()
    if len(df_cov) > 0:
        fig, ax = plt.subplots(figsize=(8, 8))
        ax.scatter(df_cov['old_mean_coverage'], df_cov['new_mean_coverage'],
                   alpha=0.5, s=20)
        min_v = min(df_cov['old_mean_coverage'].min(), df_cov['new_mean_coverage'].min()) * 0.8
        max_v = max(df_cov['old_mean_coverage'].max(), df_cov['new_mean_coverage'].max()) * 1.2
        ax.plot([min_v, max_v], [min_v, max_v], 'r--', alpha=0.5, label='y=x')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel('Old Mean Coverage')
        ax.set_ylabel('New Mean Coverage')
        ax.set_title('Mean Coverage: Old vs New')
        ax.legend()
        fig.tight_layout()
        fig.savefig(os.path.join(plot_dir, 'mean_coverage_scatter.png'), dpi=150)
        plt.close(fig)

    # 3. Assembly length scatter
    df_len = df_asm[df_asm['old_assembly_length_unambiguous'].notna() &
                     df_asm['new_assembly_length_unambiguous'].notna()].copy()
    if len(df_len) > 0:
        fig, ax = plt.subplots(figsize=(8, 8))
        ax.scatter(df_len['old_assembly_length_unambiguous'],
                   df_len['new_assembly_length_unambiguous'],
                   alpha=0.5, s=20)
        min_v = min(df_len['old_assembly_length_unambiguous'].min(),
                    df_len['new_assembly_length_unambiguous'].min()) * 0.95
        max_v = max(df_len['old_assembly_length_unambiguous'].max(),
                    df_len['new_assembly_length_unambiguous'].max()) * 1.05
        ax.plot([min_v, max_v], [min_v, max_v], 'r--', alpha=0.5, label='y=x')
        ax.set_xlabel('Old Unambiguous Length')
        ax.set_ylabel('New Unambiguous Length')
        ax.set_title('Assembly Length (Unambiguous): Old vs New')
        ax.legend()
        fig.tight_layout()
        fig.savefig(os.path.join(plot_dir, 'assembly_length_scatter.png'), dpi=150)
        plt.close(fig)

    # 4. Delta pct_ref_covered histogram
    deltas = df_asm['delta_percent_reference_covered'].dropna() * 100
    if len(deltas) > 0:
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.hist(deltas, bins=50, edgecolor='black', alpha=0.7)
        ax.axvline(0, color='r', linestyle='--', alpha=0.5)
        ax.set_xlabel('Delta % Reference Covered (New - Old)')
        ax.set_ylabel('Count')
        ax.set_title(f'Distribution of % Reference Covered Changes (n={len(deltas)})')
        fig.tight_layout()
        fig.savefig(os.path.join(plot_dir, 'delta_pct_ref_covered_hist.png'), dpi=150)
        plt.close(fig)

    # 5. Alignment identity histogram
    df_aln = df_asm[df_asm['alignment_identity'].notna()].copy()
    if len(df_aln) > 0:
        fig, ax = plt.subplots(figsize=(8, 5))
        identities = df_aln['alignment_identity'] * 100
        ax.hist(identities, bins=50, edgecolor='black', alpha=0.7)
        ax.set_xlabel('Pairwise Identity (%)')
        ax.set_ylabel('Count')
        ax.set_title(f'Assembly Identity Distribution (n={len(identities)})')
        ax.axvline(100, color='r', linestyle='--', alpha=0.3)
        fig.tight_layout()
        fig.savefig(os.path.join(plot_dir, 'alignment_identity_hist.png'), dpi=150)
        plt.close(fig)

    # 6. SNP and indel count histograms
    if len(df_aln) > 0:
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))

        snps = df_aln['snp_count'].dropna()
        if len(snps) > 0:
            max_snp = int(snps.max())
            bins_snp = range(0, max(max_snp + 2, 3))
            axes[0].hist(snps, bins=bins_snp, edgecolor='black', alpha=0.7)
            axes[0].set_xlabel('SNP Count')
            axes[0].set_ylabel('Number of Assemblies')
            axes[0].set_title(f'SNPs per Assembly (n={len(snps)})')

        indels = df_aln['indel_count_events'].dropna()
        if len(indels) > 0:
            max_indel = int(indels.max())
            bins_indel = range(0, max(max_indel + 2, 3))
            axes[1].hist(indels, bins=bins_indel, edgecolor='black', alpha=0.7)
            axes[1].set_xlabel('Indel Count (internal)')
            axes[1].set_ylabel('Number of Assemblies')
            axes[1].set_title(f'Internal Indels per Assembly (n={len(indels)})')

        fig.tight_layout()
        fig.savefig(os.path.join(plot_dir, 'snp_indel_counts.png'), dpi=150)
        plt.close(fig)

    # 7. Terminal extensions histogram
    if len(df_aln) > 0:
        ext_old = df_aln['terminal_extensions_old'].dropna()
        ext_new = df_aln['terminal_extensions_new'].dropna()
        has_ext = (ext_old > 0) | (ext_new > 0)
        if has_ext.any():
            fig, ax = plt.subplots(figsize=(8, 5))
            ax.hist(ext_new[has_ext], bins=30, alpha=0.6, label='New extends beyond Old', edgecolor='black')
            ax.hist(ext_old[has_ext], bins=30, alpha=0.6, label='Old extends beyond New', edgecolor='black')
            ax.set_xlabel('Terminal Extension (bp)')
            ax.set_ylabel('Count')
            ax.set_title('Terminal Extensions')
            ax.legend()
            fig.tight_layout()
            fig.savefig(os.path.join(plot_dir, 'terminal_extensions_hist.png'), dpi=150)
            plt.close(fig)

    # 8. Coverage vs identity
    if len(df_aln) > 0:
        df_ci = df_aln[df_aln['old_mean_coverage'].notna() & (df_aln['old_mean_coverage'] > 0)].copy()
        if len(df_ci) > 0:
            fig, ax = plt.subplots(figsize=(8, 5))
            ax.scatter(df_ci['old_mean_coverage'], df_ci['alignment_identity'] * 100,
                       alpha=0.5, s=20)
            ax.set_xscale('log')
            ax.set_xlabel('Mean Coverage')
            ax.set_ylabel('Pairwise Identity (%)')
            ax.set_title('Coverage vs Assembly Identity')
            fig.tight_layout()
            fig.savefig(os.path.join(plot_dir, 'coverage_vs_identity.png'), dpi=150)
            plt.close(fig)

    log.info(f"Generated plots in {plot_dir}")


def generate_markdown_report(df, sample_df, workspace_name, report_dir, plot_dir):
    """Generate markdown report."""
    pd, _ = get_deps()

    total_samples = len(sample_df)
    if sample_df.empty or 'old_assembly_count' not in sample_df.columns:
        samples_with_assemblies = 0
        samples_count_match = 0
    else:
        samples_with_assemblies = len(sample_df[sample_df['old_assembly_count'] > 0])
        samples_count_match = len(sample_df[sample_df['assembly_count_match']])
    samples_count_mismatch = total_samples - samples_count_match

    total_assemblies = len(df)
    if df.empty or 'alignment_identity' not in df.columns:
        df_aln = pd.DataFrame()
    else:
        df_aln = df[df['alignment_identity'].notna()]

    identical = len(df_aln[df_aln['alignment_identity'] >= 1.0]) if len(df_aln) > 0 else 0
    near_identical = len(df_aln[(df_aln['alignment_identity'] >= 0.999) & (df_aln['alignment_identity'] < 1.0)]) if len(df_aln) > 0 else 0
    minor_diff = len(df_aln[(df_aln['alignment_identity'] >= 0.99) & (df_aln['alignment_identity'] < 0.999)]) if len(df_aln) > 0 else 0
    significant_diff = len(df_aln[df_aln['alignment_identity'] < 0.99]) if len(df_aln) > 0 else 0

    if len(df_aln) > 0 and 'snp_count' in df_aln.columns:
        with_snps = len(df_aln[df_aln['snp_count'] > 0])
        with_indels = len(df_aln[df_aln['indel_count_events'] > 0])
        with_ambig = len(df_aln[df_aln['ambiguity_diffs'] > 0])
        with_terminal = len(df_aln[(df_aln['terminal_extensions_old'] > 0) | (df_aln['terminal_extensions_new'] > 0)])
        total_snps = int(df_aln['snp_count'].sum())
        total_indel_bp = int(df_aln['indel_count_bp'].sum())
        total_indel_events = int(df_aln['indel_count_events'].sum())
        total_ambig = int(df_aln['ambiguity_diffs'].sum())
        total_terminal_bp_old = int(df_aln['terminal_extensions_old'].sum())
        total_terminal_bp_new = int(df_aln['terminal_extensions_new'].sum())
        total_terminal_events_old = int(df_aln['terminal_extension_events_old'].sum())
        total_terminal_events_new = int(df_aln['terminal_extension_events_new'].sum())
    else:
        with_snps = with_indels = with_ambig = with_terminal = 0
        total_snps = total_indel_bp = total_indel_events = total_ambig = 0
        total_terminal_bp_old = total_terminal_bp_new = 0
        total_terminal_events_old = total_terminal_events_new = 0

    report_path = os.path.join(report_dir, f'report_{workspace_name}.md')
    with open(report_path, 'w') as f:
        f.write(f"# Regression Report: {workspace_name}\n\n")
        f.write(f"## Summary\n\n")
        f.write(f"| Metric | Value |\n")
        f.write(f"|--------|-------|\n")
        f.write(f"| Total samples compared | {total_samples} |\n")
        f.write(f"| Samples with assemblies (old) | {samples_with_assemblies} |\n")
        f.write(f"| Samples with matching assembly count | {samples_count_match} |\n")
        f.write(f"| Samples with mismatched assembly count | {samples_count_mismatch} |\n")
        f.write(f"| Total assembly comparisons | {total_assemblies} |\n")
        f.write(f"| Assemblies aligned | {len(df_aln)} |\n\n")

        f.write(f"## Assembly Identity\n\n")
        f.write(f"| Category | Count | % |\n")
        f.write(f"|----------|-------|---|\n")
        if len(df_aln) > 0:
            f.write(f"| Identical (100%) | {identical} | {100*identical/len(df_aln):.1f}% |\n")
            f.write(f"| Near-identical (99.9-100%) | {near_identical} | {100*near_identical/len(df_aln):.1f}% |\n")
            f.write(f"| Minor differences (99-99.9%) | {minor_diff} | {100*minor_diff/len(df_aln):.1f}% |\n")
            f.write(f"| Significant differences (<99%) | {significant_diff} | {100*significant_diff/len(df_aln):.1f}% |\n")
        f.write(f"\n")

        f.write(f"## Variant Counts\n\n")
        f.write(f"| Metric | Assemblies affected | Events | Bases |\n")
        f.write(f"|--------|--------------------:|-------:|------:|\n")
        f.write(f"| SNPs (A/C/G/T ↔ A/C/G/T) | {with_snps} | {total_snps} | {total_snps} |\n")
        f.write(f"| Internal indels | {with_indels} | {total_indel_events} | {total_indel_bp} |\n")
        f.write(f"| Ambiguity diffs (N ↔ A/C/G/T) | {with_ambig} | {total_ambig} | {total_ambig} |\n")
        f.write(f"| Terminal extensions (old only) | {with_terminal} | {total_terminal_events_old} | {total_terminal_bp_old} |\n")
        f.write(f"| Terminal extensions (new only) | {with_terminal} | {total_terminal_events_new} | {total_terminal_bp_new} |\n\n")

        # Metrics summary
        if len(df_aln) > 0:
            f.write(f"## Metrics Summary\n\n")
            f.write(f"| Metric | Median delta | Mean delta | Min | Max |\n")
            f.write(f"|--------|-------------|------------|-----|-----|\n")
            for col, label in [
                ('delta_percent_reference_covered', '% Ref Covered'),
                ('delta_mean_coverage', 'Mean Coverage'),
                ('delta_assembly_length_unambiguous', 'Unambig Length'),
            ]:
                vals = df[col].dropna()
                if len(vals) > 0:
                    f.write(f"| {label} | {vals.median():.4f} | {vals.mean():.4f} | {vals.min():.4f} | {vals.max():.4f} |\n")
            f.write(f"\n")

        # Divergent assemblies
        divergent = df_aln[df_aln['alignment_identity'] < 0.999].sort_values('alignment_identity')
        if len(divergent) > 0:
            f.write(f"## Divergent Assemblies (identity < 99.9%)\n\n")
            f.write(f"| Assembly ID | Tax Name | Identity | SNPs | Indel events (bp) | Ambig Diffs | Term Ext Old events (bp) | Term Ext New events (bp) |\n")
            f.write(f"|-------------|----------|----------|------|-------------------|-------------|--------------------------|---------------------------|\n")
            for _, row in divergent.iterrows():
                indel_ev = int(row.get('indel_count_events', 0))
                indel_bp = int(row.get('indel_count_bp', 0))
                term_old_ev = int(row.get('terminal_extension_events_old', 0))
                term_old_bp = int(row.get('terminal_extensions_old', 0))
                term_new_ev = int(row.get('terminal_extension_events_new', 0))
                term_new_bp = int(row.get('terminal_extensions_new', 0))
                f.write(f"| {row['assembly_id']} | {row['tax_name']} | "
                        f"{row['alignment_identity']*100:.3f}% | {row['snp_count']:.0f} | "
                        f"{indel_ev} ({indel_bp}) | {row['ambiguity_diffs']:.0f} | "
                        f"{term_old_ev} ({term_old_bp}) | "
                        f"{term_new_ev} ({term_new_bp}) |\n")
            f.write(f"\n")

        # Assembly count mismatches
        if 'assembly_count_match' not in sample_df.columns:
            mismatches = pd.DataFrame()
        else:
            mismatches = sample_df[~sample_df['assembly_count_match']]
        if len(mismatches) > 0:
            f.write(f"## Assembly Count Mismatches\n\n")
            f.write(f"| Sample | Old Count | New Count | Only Old | Only New |\n")
            f.write(f"|--------|-----------|-----------|----------|----------|\n")
            for _, row in mismatches.iterrows():
                f.write(f"| {row['sample_id']} | {row['old_assembly_count']} | "
                        f"{row['new_assembly_count']} | {row['assemblies_only_in_old']} | "
                        f"{row['assemblies_only_in_new']} |\n")
            f.write(f"\n")

        # Plot references
        f.write(f"## Plots\n\n")
        plot_files = sorted(os.listdir(plot_dir)) if os.path.isdir(plot_dir) else []
        for pf in plot_files:
            if pf.endswith('.png'):
                rel_plot_dir = os.path.basename(plot_dir)
                f.write(f"![{pf}]({rel_plot_dir}/{pf})\n\n")

    log.info(f"Report written to {report_path}")
    return report_path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results-dir', required=True, help='Directory with per-sample JSON results')
    parser.add_argument('--report-dir', required=True, help='Output directory for report')
    parser.add_argument('--workspace-name', required=True, help='Workspace name for report title')
    args = parser.parse_args()

    pd, plt = get_deps()

    results = load_results(args.results_dir)
    log.info(f"Loaded {len(results)} sample results")

    # Build comparison table
    comp_rows = build_comparison_table(results, args.workspace_name)
    df = pd.DataFrame(comp_rows)
    log.info(f"Total assembly comparisons: {len(df)}")

    # Build sample summary
    sample_rows = build_sample_summary(results)
    sample_df = pd.DataFrame(sample_rows)

    # Save summary TSV
    os.makedirs(args.report_dir, exist_ok=True)
    tsv_path = os.path.join(args.report_dir, f'summary_{args.workspace_name}.tsv')
    if len(df) > 0:
        df.to_csv(tsv_path, sep='\t', index=False)
        log.info(f"Summary TSV written to {tsv_path}")

    # Generate plots
    plot_dir = os.path.join(args.report_dir, f'plots_{args.workspace_name}')
    if len(df) > 0:
        generate_plots(df, plot_dir)

    # Generate markdown report
    generate_markdown_report(df, sample_df, args.workspace_name, args.report_dir, plot_dir)


if __name__ == '__main__':
    main()
