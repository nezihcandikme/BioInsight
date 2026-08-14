# Security Policy

DEConcord is a small, actively-developed research/educational tool. It is
not intended for use on sensitive or clinical data, and does not handle
credentials or secrets itself — the one exception is the optional AI
summary feature (`explain_results=True`), which sends already-computed
differential expression results to the Anthropic API using an
`ANTHROPIC_API_KEY` you provide via your own environment. DEConcord never
stores, logs, or transmits that key anywhere else.

## Supported versions

This project is pre-1.0 and moves quickly. Only the latest released
version is supported; please update before reporting an issue.

## Reporting a vulnerability

If you find a security issue (e.g. unsafe deserialization of untrusted
input, a path traversal in file loading, or a similar concrete
vulnerability — not a general bug), please email
[nezihcandikme@gmail.com](mailto:nezihcandikme@gmail.com) directly rather
than opening a public issue. Include a minimal reproduction if possible.
Given this is a one-person project, response times won't be
enterprise-grade, but reports will be taken seriously and credited.

For ordinary bugs (incorrect statistics, crashes on valid input,
documentation errors), please use the normal
[issue tracker](https://github.com/nezihcandikme/BioInsight/issues)
instead.
