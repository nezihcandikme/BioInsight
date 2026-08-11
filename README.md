BioInsight

BioInsight is a Python pipeline for taking an RNA-seq count matrix from "okay, I have a CSV" to actual exploratory results: validation, sample quality control, normalization, differential expression, plots, and pathway enrichment — with an optional layer that asks an LLM to explain the output in plain language.

I'm building it because I wanted to understand what actually happens between getting biological data and claiming that it means something. Turns out there are approximately seventeen ways to produce convincing nonsense before breakfast, so BioInsight is currently focused on making those mistakes loud, testable, and difficult to ignore.

Current version: v0.4.0. It works, it has tests, and it is actively improving. It is also an educational project — not a production replacement for DESeq2, edgeR, or an actual bioinformatician who has seen your experimental design.

So what is this, exactly?

BioInsight started with one extremely glamorous problem: validating a CSV.

Then the CSV raised questions.

Are these really non-negative integer counts?

Are genes rows and samples columns, or did someone rotate reality by 90 degrees?

Does one sample have dramatically fewer reads than the others?

How do I compare samples with different sequencing depths?

Which genes differ between conditions?

If I test thousands of genes, how many "discoveries" are just statistical noise wearing a lab coat?

Do those genes converge on any interesting pathways?

Each question became a piece of the pipeline. The result follows the same broad reasoning order a real analysis does: load → validate → inspect → normalize → compare → visualize → interpret.

The code running is step one. The numbers meaning what I think they mean is the actual objective.

Why I'm building it

I'm a high-school student learning programming, statistics, and computational biology by building things that are slightly beyond what I currently know how to build.

Before BioInsight, I understood RNA-seq in the dangerously comfortable summary-version: count reads per gene, compare two conditions, find changed genes. That explanation is technically related to reality, but it skips nearly everything capable of ruining an analysis.

Building the pipeline forced me to confront those missing layers directly: what makes a count matrix valid, why library size matters, why normalization is not one universal operation, why a difference of means only becomes a fold change on the right scale, why thousands of simultaneous tests require correction, and why statistical significance is not the same thing as biological importance.

BioInsight is the record of learning those lessons in code — one function, one broken test, and one increasingly specific error message at a time. The full decision-by-decision version of that record lives in `DEVLOG.md`.

What it does, roughly

Raw counts go in. Along the way: the data gets checked for the kind of problems that quietly wreck an analysis, each sample gets sanity-checked against the others, samples get put on a comparable scale, genes get tested for meaningful differences between conditions, the results get corrected for the fact that testing thousands of things at once produces false alarms, and — if the numbers hold up — the changed genes get checked for whether they cluster into any known biological pathways. Plots come out along the way so you can actually look at the data instead of just trusting a table.

Every piece of that also works on its own, in case you want to stop and poke at an intermediate result instead of running the whole thing end to end.

Where it stands

It works, on real-shaped toy data, with a test for every rough edge I've found so far (and I keep finding more). What it isn't yet: validated against the field's actual gold-standard tools, so treat its output as a fast first look, not a final answer. It also only understands one kind of data right now — bulk RNA-seq count matrices. Broader ambitions exist; pretending they're already real would just make the README more advanced than the software.

Where it's headed

The long-term idea is bigger than RNA-seq: something that can look at a scientific dataset, figure out what kind of analysis actually applies to it, ask for whatever context is missing, run validated methods instead of guessing, and explain the result honestly instead of impressively.

That's a long way off from "upload arbitrary CSV, receive truth," so the near-term plan stays concrete: check the current results against DESeq2 and edgeR on a real public dataset, make pathway enrichment accept standard gene-set formats, and only then start adding the more ambitious layers. Barely-expressed genes no longer get tested by default reasoning alone — that filter exists now, you just have to turn it on.

Current mission, in short: make the statistical layer harder to fool before teaching it new tricks.

The detailed version — what changed, when, and why — is in `DEVLOG.md`.
