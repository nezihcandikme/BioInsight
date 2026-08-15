import pandas as pd
import pytest

from deconcord.concordance.threshold_sensitivity import compute_threshold_sensitivity


def _results_a():
    return pd.DataFrame({
        # Varied (not constant) so Pearson/Spearman on log fold change
        # is well-defined -- irrelevant to what this module reports, but
        # avoids a scipy ConstantInputWarning inside compute_de_concordance.
        "log_fold_change": [1.0, 1.2, 0.9],
        # g1: significant at every alpha tried. g2: significant only at
        # the looser thresholds (0.05, 0.1), not at 0.01. g3: never
        # significant.
        "adjusted_p_value": [0.005, 0.03, 0.5],
    }, index=["g1", "g2", "g3"])


def _results_b():
    return pd.DataFrame({
        "log_fold_change": [1.1, 0.8, 1.3],
        # g1: significant at every alpha, same as A (stable and
        # concordant). g2: never significant in B (A's borderline call
        # doesn't even show up here). g3: significant at every alpha,
        # the opposite of A (stable in each table on its own, but the
        # two tables disagree at every threshold).
        "adjusted_p_value": [0.005, 0.5, 0.005],
    }, index=["g1", "g2", "g3"])


def test_compute_threshold_sensitivity_by_alpha_table():
    result = compute_threshold_sensitivity(
        _results_a(), _results_b(), alphas=[0.01, 0.05, 0.1], name_a="A", name_b="B",
    )
    by_alpha = result["by_alpha"]

    assert list(by_alpha.index) == [0.01, 0.05, 0.1]
    assert list(by_alpha["genes_compared"]) == [3, 3, 3]
    # alpha=0.01: sig_A={g1}, sig_B={g1,g3} -> jaccard = 1/2
    assert by_alpha.loc[0.01, "jaccard_index"] == pytest.approx(1 / 2)
    # alpha=0.05 and 0.1: sig_A={g1,g2}, sig_B={g1,g3} -> jaccard = 1/3
    assert by_alpha.loc[0.05, "jaccard_index"] == pytest.approx(1 / 3)
    assert by_alpha.loc[0.1, "jaccard_index"] == pytest.approx(1 / 3)
    assert list(by_alpha["significant_in_both"]) == [1, 1, 1]


def test_compute_threshold_sensitivity_gene_stability_fractions():
    result = compute_threshold_sensitivity(
        _results_a(), _results_b(), alphas=[0.01, 0.05, 0.1], name_a="A", name_b="B",
    )
    stability = result["gene_stability"]

    assert stability.loc["g1", "frac_significant_A"] == pytest.approx(1.0)
    assert stability.loc["g1", "frac_significant_B"] == pytest.approx(1.0)
    assert stability.loc["g2", "frac_significant_A"] == pytest.approx(2 / 3)
    assert stability.loc["g2", "frac_significant_B"] == pytest.approx(0.0)
    assert stability.loc["g3", "frac_significant_A"] == pytest.approx(0.0)
    assert stability.loc["g3", "frac_significant_B"] == pytest.approx(1.0)


def test_compute_threshold_sensitivity_stable_vs_sensitive_genes():
    result = compute_threshold_sensitivity(
        _results_a(), _results_b(), alphas=[0.01, 0.05, 0.1], name_a="A", name_b="B",
    )

    # g1 and g3 both have a fraction of 0.0 or 1.0 in *each* table
    # individually, so the significance call itself never flips with
    # threshold -- g3 counts as "stable" here even though A and B
    # disagree with each other about g3 at every alpha. Stability is
    # about surviving a threshold change, not about the two tables
    # agreeing with each other.
    assert list(result["stable_genes"]) == ["g1", "g3"]
    # g2 flips between significant and not significant somewhere in the
    # swept range in table A, so its concordance status with B depends
    # on which threshold you picked.
    assert list(result["sensitive_genes"]) == ["g2"]


def test_compute_threshold_sensitivity_single_alpha_matches_compute_de_concordance():
    from deconcord.concordance.methods import compute_de_concordance

    a, b = _results_a(), _results_b()
    direct = compute_de_concordance(a, b, name_a="A", name_b="B", alpha=0.05)
    swept = compute_threshold_sensitivity(a, b, alphas=[0.05], name_a="A", name_b="B")

    assert swept["by_alpha"].loc[0.05, "jaccard_index"] == pytest.approx(
        direct["summary"]["jaccard_index"]
    )
    assert swept["by_alpha"].loc[0.05, "significant_in_both"] == direct["summary"]["significant_in_both"]


def test_compute_threshold_sensitivity_empty_alphas_raises():
    with pytest.raises(ValueError, match="alphas must not be empty"):
        compute_threshold_sensitivity(_results_a(), _results_b(), alphas=[])


def test_compute_threshold_sensitivity_duplicate_alphas_raises():
    with pytest.raises(ValueError, match="duplicates"):
        compute_threshold_sensitivity(_results_a(), _results_b(), alphas=[0.05, 0.05])


def test_compute_threshold_sensitivity_invalid_alpha_raises():
    with pytest.raises(ValueError, match=r"\(0, 1\]"):
        compute_threshold_sensitivity(_results_a(), _results_b(), alphas=[0.05, 0])

    with pytest.raises(ValueError, match=r"\(0, 1\]"):
        compute_threshold_sensitivity(_results_a(), _results_b(), alphas=[0.05, 1.5])


def test_compute_threshold_sensitivity_propagates_compute_de_concordance_errors():
    a = pd.DataFrame({"log_fold_change": [1.0]}, index=["g1"])  # no p-value column
    b = _results_b()

    with pytest.raises(ValueError, match="missing column"):
        compute_threshold_sensitivity(a, b, name_a="A")


def test_compute_threshold_sensitivity_alphas_order_independent():
    # Passing alphas out of order shouldn't change the result -- the
    # function sorts them internally before sweeping.
    a, b = _results_a(), _results_b()
    ascending = compute_threshold_sensitivity(a, b, alphas=[0.01, 0.05, 0.1], name_a="A", name_b="B")
    shuffled = compute_threshold_sensitivity(a, b, alphas=[0.1, 0.01, 0.05], name_a="A", name_b="B")

    pd.testing.assert_frame_equal(ascending["by_alpha"], shuffled["by_alpha"])
    pd.testing.assert_frame_equal(ascending["gene_stability"], shuffled["gene_stability"])
