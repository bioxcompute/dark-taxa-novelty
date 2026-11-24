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
| `baseline.py` | the baseline below, in one file, numpy and nothing else |

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

Counting four-letter words. Take the 4-mer profile of a query, compare it to one
centroid per reference group, and keep the nearest. No training, no weights, no
GPU, and it runs in about twenty seconds on a laptop.

`baseline.py` in this repository is that method, in one file, with numpy as its
only dependency. Every number below came out of it, measured on the files beside
it on 24 November 2026, so you can reproduce the row you want to beat before you
try to beat it.

| method | length | AUC | 95% over rows | 95% over taxa |
|---|---|---|---|---|
| one centroid per group | 300 bases | **0.7329** | [0.6976, 0.7674] | [0.6249, 0.8369] |
| one centroid per group | 512 bases | **0.7578** | [0.7225, 0.7910] | [0.6426, 0.8680] |
| one centroid per group | full read | **0.7741** | [0.7389, 0.8075] | [0.6578, 0.8874] |
| nearest single sequence | full read | 0.8022 | [0.7690, 0.8350] | [0.6727, 0.9174] |

Guessing is 0.5000. The median sequence here is about 1,450 bases, so "full
read" is roughly three times the longest crop.

**More context helps, and it helps monotonically**: 0.7329, 0.7578, 0.7741 as
the crop goes from 300 to 512 to the whole read, with nothing else changed. That
is why a model read at one length may only be compared with a baseline measured
at that same length. An extra 212 bases is worth about 0.025 AUC here, which is
larger than the gap between several published models.

**One centroid per group, not the nearest single reference sequence.** Taking the
nearest single sequence rewards a query that happens to resemble one odd record,
and gives a large group more chances to win simply by being large. Measured, that
choice is worth 0.0281 and it goes *up*, not down: 0.8022 against 0.7741. So the
row we ask you to beat is the more conservative of the two, and `--nearest` runs
it the other way if you want to see that for yourself.

An earlier version of this file said 0.755 for this baseline and put the centroid
choice at "about 0.09". Neither survived being measured on the published files.
The number was close to the 512-base crop and the 0.09 was three times too large,
so this paragraph is here instead of a quiet edit.

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
