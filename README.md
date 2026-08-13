OmicForge

OmicForge is a Python pipeline for taking an RNA-seq count matrix from "okay, I have a CSV" to actual exploratory results: validation, sample quality control, normalization, differential expression, plots, and pathway enrichment — with an optional layer that asks an LLM to explain the output in plain language.

I'm building it because I wanted to understand what actually happens between getting biological data and claiming that it means something. Turns out there are approximately seventeen ways to produce convincing nonsense before breakfast, so OmicForge is currently focused on making those mistakes loud, testable, and difficult to ignore.

Current version: v0.9.0. It works, it has tests, and it is actively improving. It is also an educational project — not a production replacement for DESeq2, edgeR, or an actual bioinformatician who has seen your experimental design.

So what is this, exactly?

OmicForge started with one extremely glamorous problem: validating a CSV.

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

Before OmicForge, I understood RNA-seq in the dangerously comfortable summary-version: count reads per gene, compare two conditions, find changed genes. That explanation is technically related to reality, but it skips nearly everything capable of ruining an analysis.

Building the pipeline forced me to confront those missing layers directly: what makes a count matrix valid, why library size matters, why normalization is not one universal operation, why a difference of means only becomes a fold change on the right scale, why thousands of simultaneous tests require correction, and why statistical significance is not the same thing as biological importance.

OmicForge is the record of learning those lessons in code — one function, one broken test, and one increasingly specific error message at a time. The full decision-by-decision version of that record lives in `DEVLOG.md`.

What it does, roughly

Raw counts go in. Along the way: the data gets checked for the kind of problems that quietly wreck an analysis, each sample gets sanity-checked against the others, samples get put on a comparable scale, genes get tested for meaningful differences between conditions (either gene-by-gene, or with a second method that borrows statistical power across genes when a gene's own replicates are too noisy to trust alone), the results get corrected for the fact that testing thousands of things at once produces false alarms, and — if the numbers hold up — the changed genes get checked for whether they cluster into any known biological pathways. Plots come out along the way so you can actually look at the data instead of just trusting a table.

Every piece of that also works on its own, in case you want to stop and poke at an intermediate result instead of running the whole thing end to end. It also runs as a command-line tool now, not just as a library you import — point it at a CSV and two group names and it writes the results table, the plots, and the enrichment table to a folder, without needing a Python script in between. Results can also come back with actual gene names attached instead of just the accession IDs they were tested under, if you hand it a mapping — nothing fetched over the network, just a local file you already have. Pathway enrichment now has a second, opt-in source too: alongside the original local GMT-file lookup, it can query g:Profiler's live, curated databases directly, if you're running somewhere with real internet access. And if the data you want to analyze lives in a GEO series rather than on your disk already, there's now a small set of helpers that fetch a series' supplementary files and sample metadata straight from NCBI, so getting from a GSE accession to a usable count matrix doesn't mean hand-building the FTP URL yourself.

Where it stands

It works, on real-shaped toy data and on a real public dataset, with a test for every rough edge I've found so far (and I keep finding more). It's also now actually been checked against DESeq2 and edgeR, not just labeled "unvalidated" and left there: on a real RNA-seq dataset, its fold-change estimates track both tools closely, and every gene it calls significant, both of the field-standard tools also call significant — no false alarms, with either testing method. What the default, simplest method doesn't do is find most of what DESeq2 and edgeR find; a plain per-gene t-test has measurably less statistical power than models built to borrow strength across genes. The second, opt-in method exists because of that exact finding — it borrows that same kind of strength across genes, and on the same real dataset it found more than double the significant genes the default method did, while keeping the same zero-false-alarm result. Neither number is the final word; both are in `benchmarks/`, along with the methodology and the honest caveats. It also only understands one kind of data right now — bulk RNA-seq count matrices. Broader ambitions exist; pretending they're already real would just make the README more advanced than the software.

Where it's headed

The long-term idea is bigger than RNA-seq: something that can look at a scientific dataset, figure out what kind of analysis actually applies to it, ask for whatever context is missing, run validated methods instead of guessing, and explain the result honestly instead of impressively.

That's a long way off from "upload arbitrary CSV, receive truth." The concrete near-term items keep clearing: barely-expressed genes no longer get tested by default reasoning alone, pathway enrichment reads the same .gmt files MSigDB ships, and the results are now checked against DESeq2 and edgeR instead of just described as unchecked — all three used to be roadmap items, none of them still are.

Current mission, in short: keep finding out exactly where the statistical layer's numbers can and can't be trusted, before teaching it new tricks. The moderated testing method is the first real example of that loop actually closing — a measured weakness turned into a specific fix, then checked against the same benchmark that found the weakness in the first place.

The detailed version — what changed, when, and why — is in `DEVLOG.md`.
