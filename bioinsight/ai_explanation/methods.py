import anthropic
import os
import pandas as pd
from dotenv import load_dotenv


def explain_de_results(results_df: pd.DataFrame) -> str:
    """
    Generates an explanation of differential expression results
    using the Anthropic Claude API.

    Args:
        results_df (pd.DataFrame): A DataFrame containing differential
            expression results.

    Returns:
        str: The AI-generated explanation of the differential expression results.
    """
    load_dotenv()

    api_key = os.environ["ANTHROPIC_API_KEY"]
    client = anthropic.Anthropic(api_key=api_key)

    n_total = len(results_df)
    n_significant = results_df["significant"].sum()

    top_genes = (
        results_df
        .sort_values("adjusted_p_value")
        .head(5)
    )

    prompt = f"""
I am summarizing the results of an RNA-seq differential expression analysis.

Total number of genes tested: {n_total}
Number of statistically significant genes: {n_significant}

Top 5 most significant genes:
{top_genes[['log_fold_change', 'adjusted_p_value']].to_string()}

Briefly and clearly explain these results to a researcher.
Describe the overall pattern of the analysis and mention any important
considerations or potential issues that should be kept in mind when
interpreting the results.
"""

    response = client.messages.create(
        model="claude-sonnet-4-5",
        max_tokens=1024,
        messages=[{"role": "user", "content": prompt}],
    )

    return response.content[0].text