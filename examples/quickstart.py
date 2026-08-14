"""
DEConcord quick-start example.

Runs the whole intended workflow end to end on a small, synthetic,
deterministic dataset: load a count matrix, validate it, run differential
expression two different (reasonable) ways, and check how much the
conclusions actually agree.

No network access, no R, no external data download -- everything here
runs from what's already in this repository, and every number it prints
is reproducible run to run (fixed random seed).

Usage:
    python examples/quickstart.py

Outputs (written to examples/output/, created if missing):
    quickstart_counts.csv   the synthetic input count matrix
    volcano_welch.png       volcano plot from the Welch's t-test DE run
    concordance_summary.txt human-readable concordance report
"""

from pathlib import Path

import numpy as np
import pandas as pd

import deconcord as dc

OUTPUT_DIR = Path(__file__).parent / "output"
N_GENES = 3000
N_TRUE_DE_GENES = 60
N_SAMPLES_PER_GROUP = 4
RANDOM_SEED = 0


def make_synthetic_counts() -> pd.DataFrame:
    """
    Build a small, deterministic synthetic count matrix: most genes have
    no real difference between the two groups (noise only), and a known
    subset of genes have a real, planted log2 fold change -- so this
    script's own DE calls have a known-correct answer to be checked
    against, not just "some numbers came out."
    """
    rng = np.random.default_rng(RANDOM_SEED)
    samples = [f"sample_{i + 1}" for i in range(2 * N_SAMPLES_PER_GROUP)]
    group_1 = samples[:N_SAMPLES_PER_GROUP]
    group_2 = samples[N_SAMPLES_PER_GROUP:]

    baseline = rng.uniform(low=20, high=500, size=N_GENES)
    is_de = np.zeros(N_GENES, dtype=bool)
    is_de[:N_TRUE_DE_GENES] = True
    rng.shuffle(is_de)
    true_log2fc = np.where(is_de, rng.choice([-1, 1], size=N_GENES) * rng.uniform(1.5, 3.0, size=N_GENES), 0.0)

    counts = {}
    for sample in group_1:
        counts[sample] = rng.poisson(lam=baseline).astype(int)
    for sample in group_2:
        shifted_mean = baseline * (2.0**true_log2fc)
        counts[sample] = rng.poisson(lam=shifted_mean).astype(int)

    gene_ids = [f"gene_{i:04d}" for i in range(N_GENES)]
    return pd.DataFrame(counts, index=gene_ids), group_1, group_2


def main() -> None:
    OUTPUT_DIR.mkdir(exist_ok=True)

    print(f"DEConcord {dc.__version__} quick-start example\n")

    # 1. Generate and write a synthetic count matrix, then load + validate
    #    it the same way any real input would go through the pipeline --
    #    this is what load_count_matrix / validate_counts actually check.
    counts_df, group_1, group_2 = make_synthetic_counts()
    counts_path = OUTPUT_DIR / "quickstart_counts.csv"
    counts_df.to_csv(counts_path, index_label="gene_id")

    counts = dc.load_count_matrix(str(counts_path))
    print(f"Loaded and validated {counts.shape[0]} genes x {counts.shape[1]} samples from {counts_path.name}")

    # 2. Normalize before differential expression, as any of DEConcord's
    #    DE functions expect.
    cpm = dc.compute_cpm(counts)
    log2_cpm = dc.log2_transform(cpm)

    # 3. Run differential expression two different, reasonable ways: a
    #    plain Welch's t-test and the empirical-Bayes moderated t-test.
    #    This is exactly the kind of "reasonable analytical decision"
    #    DEConcord exists to stress-test -- same data, same groups,
    #    different (both legitimate) significance test.
    welch_results = dc.run_differential_expression(log2_cpm, group_1, group_2, method="welch")
    moderated_results = dc.run_differential_expression(log2_cpm, group_1, group_2, method="moderated")

    print(f"Welch's t-test:      {welch_results['significant'].sum()} significant genes")
    print(f"Moderated t-test:    {moderated_results['significant'].sum()} significant genes")

    # 4. Ask the actual question DEConcord is for: how much do these two
    #    reasonable choices agree with each other?
    concordance = dc.compute_de_concordance(
        welch_results, moderated_results, name_a="welch", name_b="moderated",
    )
    summary = concordance["summary"]

    report_lines = [
        f"DEConcord {dc.__version__} quick-start -- method concordance summary",
        "=" * 60,
        f"Genes compared:          {summary['genes_compared']}",
        f"Significant (welch):     {summary['significant_welch']}",
        f"Significant (moderated): {summary['significant_moderated']}",
        f"Significant in both:     {summary['significant_in_both']}",
        f"Jaccard index:           {summary['jaccard_index']:.3f}",
        f"Pearson r (log2FC):      {summary['pearson_r_lfc']:.3f}",
        f"Directional agreement:   {summary['directional_agreement']:.3f}",
        f"Discordant genes:        {len(concordance['discordant_genes'])}",
        "",
        "A high Jaccard index and directional agreement here mean the",
        "significant-gene conclusion is not an artifact of picking one",
        "particular significance test over the other. See METHODOLOGY.md",
        "for what these numbers do and don't mean.",
    ]
    report = "\n".join(report_lines)
    print("\n" + report)

    report_path = OUTPUT_DIR / "concordance_summary.txt"
    report_path.write_text(report + "\n")

    # 5. A figure, because a table of numbers is harder to sanity-check
    #    at a glance than a plot.
    fig = dc.plot_volcano(welch_results)
    fig_path = OUTPUT_DIR / "volcano_welch.png"
    fig.savefig(fig_path, dpi=150, bbox_inches="tight")

    print(f"\nWrote {counts_path}, {report_path}, and {fig_path}")


if __name__ == "__main__":
    main()
