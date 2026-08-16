"""Resampling stability: rerunning DE on perturbed sample sets and checking
how often each gene's significance call reappears. See
src/deconcord/concordance/resampling_stability.py's module docstring for
why bootstrap resampling isn't offered (only "subsample" and
"leave_one_out" are).
"""
from fractions import Fraction

import pandas as pd
import pytest

from deconcord.concordance.resampling_stability import compute_resampling_stability
from deconcord.differential_expression.methods import run_differential_expression

GROUP_1 = ["s1", "s2", "s3"]
GROUP_2 = ["s4", "s5", "s6"]


def _df():
    # "strong": huge, low-noise gap between groups -- stays significant no
    # matter which sample is dropped or subsampled out.
    # "borderline": significant with every sample present, but losing
    # either s1 or s2 (both slightly pull group_1's mean down) flips it to
    # non-significant. Losing anything else leaves it significant. Values
    # and the resulting flip pattern were checked against a real run of
    # run_differential_expression before being hardcoded here.
    # "flat": near-identical means in both groups -- never significant.
    return pd.DataFrame({
        "s1": [1.0, 5.0, 5.0],
        "s2": [1.2, 5.05, 5.0],
        "s3": [0.8, 4.5, 5.0],
        "s4": [9.0, 3.0, 5.02],
        "s5": [9.3, 3.05, 4.98],
        "s6": [8.7, 3.1, 5.03],
    }, index=["strong", "borderline", "flat"])


def test_leave_one_out_baseline_matches_direct_call():
    df = _df()
    result = compute_resampling_stability(df, GROUP_1, GROUP_2, resample_method="leave_one_out")
    direct = run_differential_expression(df, GROUP_1, GROUP_2, alpha=0.05, lfc_threshold=1.0, method="welch")
    pd.testing.assert_frame_equal(result["baseline_results"], direct)


def test_leave_one_out_gene_stability_fractions():
    result = compute_resampling_stability(_df(), GROUP_1, GROUP_2, resample_method="leave_one_out")
    stability = result["gene_stability"]

    assert stability.loc["strong", "baseline_significant"]
    assert stability.loc["strong", "frac_significant"] == pytest.approx(1.0)
    assert stability.loc["borderline", "baseline_significant"]
    assert stability.loc["borderline", "frac_significant"] == pytest.approx(float(Fraction(4, 6)))
    assert not stability.loc["flat", "baseline_significant"]
    assert stability.loc["flat", "frac_significant"] == pytest.approx(0.0)


def test_leave_one_out_stable_vs_sensitive_genes():
    result = compute_resampling_stability(_df(), GROUP_1, GROUP_2, resample_method="leave_one_out")

    # "strong" replicates every time (stable, baseline-significant),
    # "flat" is never significant either way (stable, baseline-non-
    # significant), "borderline" only replicates 4/6 of the time, below
    # the default stability_threshold of 0.9, so it's sensitive despite
    # being significant in the baseline run.
    assert list(result["stable_genes"]) == ["flat", "strong"]
    assert list(result["sensitive_genes"]) == ["borderline"]


def test_leave_one_out_summary():
    result = compute_resampling_stability(_df(), GROUP_1, GROUP_2, resample_method="leave_one_out")
    summary = result["summary"]

    assert summary["resample_method"] == "leave_one_out"
    assert summary["n_iterations"] == 6  # 3 + 3 samples, each dropped once
    assert summary["n_baseline_significant"] == 2  # strong, borderline
    assert summary["mean_jaccard_to_baseline"] == pytest.approx(float(Fraction(5, 6)))
    assert summary["baseline_replication_rate"] == pytest.approx(float(Fraction(5, 6)))


def test_leave_one_out_needs_at_least_three_samples_per_group():
    df = _df()
    with pytest.raises(ValueError, match="leave_one_out needs at least 3 samples"):
        compute_resampling_stability(df, ["s1", "s2"], GROUP_2, resample_method="leave_one_out")


def test_subsample_reproducible_with_same_random_state():
    df = _df()
    a = compute_resampling_stability(
        df, GROUP_1, GROUP_2, resample_method="subsample", n_iterations=50, random_state=7,
    )
    b = compute_resampling_stability(
        df, GROUP_1, GROUP_2, resample_method="subsample", n_iterations=50, random_state=7,
    )
    pd.testing.assert_frame_equal(a["gene_stability"], b["gene_stability"])
    assert a["summary"] == b["summary"]


def test_subsample_strong_and_flat_genes_are_stable():
    result = compute_resampling_stability(
        _df(), GROUP_1, GROUP_2, resample_method="subsample", n_iterations=100, random_state=1,
    )
    stability = result["gene_stability"]

    # subsample_fraction defaults to 0.8; with 3 samples per group that
    # rounds to keeping 2, the same "drop one" gap "strong"/"flat" were
    # designed to survive regardless of which one is dropped.
    assert stability.loc["strong", "frac_significant"] == pytest.approx(1.0)
    assert stability.loc["flat", "frac_significant"] == pytest.approx(0.0)
    assert "strong" in result["stable_genes"]
    assert "flat" in result["stable_genes"]


def test_subsample_honors_n_iterations():
    result = compute_resampling_stability(
        _df(), GROUP_1, GROUP_2, resample_method="subsample", n_iterations=37, random_state=1,
    )
    assert result["summary"]["n_iterations"] == 37
    assert result["summary"]["resample_method"] == "subsample"


def test_invalid_resample_method_raises():
    with pytest.raises(ValueError, match="resample_method must be one of"):
        compute_resampling_stability(_df(), GROUP_1, GROUP_2, resample_method="bootstrap")


def test_invalid_stability_threshold_raises():
    with pytest.raises(ValueError, match="stability_threshold must be in"):
        compute_resampling_stability(_df(), GROUP_1, GROUP_2, stability_threshold=0.5)
    with pytest.raises(ValueError, match="stability_threshold must be in"):
        compute_resampling_stability(_df(), GROUP_1, GROUP_2, stability_threshold=1.1)


def test_invalid_n_iterations_raises():
    with pytest.raises(ValueError, match="n_iterations must be a positive integer"):
        compute_resampling_stability(_df(), GROUP_1, GROUP_2, resample_method="subsample", n_iterations=0)
    with pytest.raises(ValueError, match="n_iterations must be a positive integer"):
        compute_resampling_stability(_df(), GROUP_1, GROUP_2, resample_method="subsample", n_iterations=-5)


def test_invalid_subsample_fraction_raises():
    with pytest.raises(ValueError, match="subsample_fraction must be in"):
        compute_resampling_stability(_df(), GROUP_1, GROUP_2, resample_method="subsample", subsample_fraction=0)
    with pytest.raises(ValueError, match="subsample_fraction must be in"):
        compute_resampling_stability(_df(), GROUP_1, GROUP_2, resample_method="subsample", subsample_fraction=1)


def test_subsample_fraction_leaving_too_few_samples_raises():
    df = _df()
    with pytest.raises(ValueError, match="fewer than 2 samples"):
        compute_resampling_stability(
            df, ["s1", "s2"], ["s4", "s5"], resample_method="subsample", subsample_fraction=0.4,
        )


def test_propagates_run_differential_expression_errors():
    df = _df()
    with pytest.raises(ValueError, match="not present in the count matrix"):
        compute_resampling_stability(df, ["s1", "nonexistent"], GROUP_2)
