---
license: cc-by-4.0
task_categories:
  - zero-shot-classification
tags:
  - biology
  - genomics
  - eDNA
  - 18S
  - 16S
  - benchmark
  - novelty-detection
pretty_name: Dark taxa novelty, v1.0
size_categories:
  - 1K<n<10K
---

# Dark taxa novelty, v1.0

**Can a model recognise a taxon it has never seen?**

Below five hundred metres, most of what you sequence from seawater matches no
named species. A reference database answers that with silence, because a
catalogue can only return what is in it. A model that reads sequence can in
principle say "this looks like nothing I know" without knowing what it is
instead. This set measures whether it actually can.

There was no fixed set for that question, so everyone reported it on their own
data and no two numbers could be compared. This one is frozen and versioned, so
a number measured today still means something next year.

## What is in it

| file | what |
|---|---|
| `reference.fasta` | 960 sequences, 16 from each of 60 reference groups |
| `reference.groups` | the group label for each, one per line, same order |
| `known.fasta` | 400 queries from those same groups, not in the reference set |
| `novel.fasta` | 400 queries from 40 taxa held out entirely |
| `benchmark.json` | the held-out taxa, the counts, and the exact split |

Every sequence carries its accession, its source database and its taxon in the
FASTA header. 1,665 of 1,760 could be traced back to an accession; the rest
carry a positional identifier instead.

## How to use it

Score each query with your model, in whatever way gives one number per sequence
where higher means more expected. Then compute the AUC: how often a sequence
from `known.fasta` outranks one from `novel.fasta`. 0.5 is guessing.

Report a bootstrap interval with it. An AUC without one invites reading noise as
progress.

Two things to keep the comparison fair. Give every model the same sequences and
the same length, because more context is easier and that advantage has nothing
to do with the model. And compare rankings, not likelihoods: the log-likelihoods
of two models sit on different scales.

## The baseline to beat

Counting four-letter words. Take the 4-mer profile of a query and compare it to
one centroid per reference group, then take the nearest.

**One centroid per group, not the nearest single reference sequence.** That
choice moves the number by about 0.09, which is more than most models differ by.
Taking the nearest single sequence rewards a query that happens to resemble one
odd record, and gives a large group more chances to win simply by being large.

On this set, that baseline reaches **AUC 0.755, interval [0.720, 0.786]**.

If you see a different number for counting, the difference is in how the
sequences are cropped, not in the model. We report 0.774 elsewhere for the same
baseline on full-length reads; that is the same method on longer input. Crop
every model the same way and the comparison holds.

## How the split was made

Whole taxa are held out, never random rows. Sequences within one taxon are not
independent of each other, so a random split lets a model answer from memory and
call it recognition.

The forty held-out taxa are listed in `benchmark.json`. They were chosen by seed
from the groups with between 30 and 4,000 records, so every one has enough
sequences to measure with and none is so large that holding it out reshapes the
rest.

## Sources and licence

Sequences come from:

- [SILVA](https://www.arb-silva.de), CC-BY 4.0
- [PR2](https://pr2-database.org), MIT

Both permit redistribution and commercial use with attribution, which is given
here and in the FASTA headers. This set is released under CC-BY 4.0 on the same
terms.

## Versioning

A benchmark that changes is not a benchmark. This set is frozen at v1.0. If a
second version is ever needed it will sit beside this one, not replace it.

---

Built by [BioX](https://bioxcompute.com). Code:
[github.com/bioxcompute/eh14x](https://github.com/bioxcompute/eh14x)
