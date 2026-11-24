"""The baseline to beat: counting four-letter words.

No training, no weights, no GPU. It builds one 4-mer profile per reference
group, scores a query by its cosine similarity to the nearest group, and asks
how often a sequence from a group the set knows outranks one from a taxon held
out entirely.

    python baseline.py                 # on the sequences as they are
    python baseline.py --crop 300      # every sequence cut to 300 bases

WHY ONE CENTROID PER GROUP AND NOT THE NEAREST SINGLE SEQUENCE
Taking the nearest single reference sequence rewards a query that happens to
resemble one odd record, and it gives a large group more chances to win simply
by being large. That choice moves the number more than most models differ by,
so it is stated rather than buried. `--nearest` runs it the other way so you
can see the size of it yourself.

WHY A BOOTSTRAP INTERVAL, AND WHY OVER TAXA
An AUC without an interval invites reading noise as progress. Resampling rows
treats two sequences from one taxon as two independent observations, and they
are not: whole taxa were held out, so the taxon is the unit. Both intervals are
printed, because the narrower one is what most papers show.

Needs numpy and nothing else.
"""
from __future__ import annotations

import argparse
import itertools
import sys

import numpy as np

K = 4
MERS = {"".join(p): i for i, p in enumerate(itertools.product("ACGT", repeat=K))}


def lees(pad: str) -> list[tuple[str, str]]:
    uit: list[list[str]] = []
    for regel in open(pad):
        if regel.startswith(">"):
            uit.append([regel[1:].strip(), ""])
        elif uit:
            uit[-1][1] += regel.strip()
    return [(k, s.upper()) for k, s in uit]


def profiel(seq: str) -> np.ndarray:
    """Hoe vaak elk van de 256 vierletterwoorden voorkomt, genormaliseerd.

    Een woord met een N of een andere letter erin wordt overgeslagen en niet
    naar A geraden: dat zou een telling verzinnen die er niet is.
    """
    v = np.zeros(len(MERS), dtype=np.float64)
    for i in range(len(seq) - K + 1):
        j = MERS.get(seq[i:i + K])
        if j is not None:
            v[j] += 1.0
    n = np.linalg.norm(v)
    return v / n if n else v


def auc(hoog: np.ndarray, laag: np.ndarray) -> float:
    """Hoe vaak een score uit `hoog` boven een uit `laag` staat, gelijkspel half.

    Met de rangen gerekend en niet met alle paren langs elkaar, zodat het ook
    bij grote sets nog in een oogwenk klaar is. Dit is dezelfde maat als de
    oppervlakte onder de ROC-curve.
    """
    alles = np.concatenate([hoog, laag])
    r = np.empty(len(alles), dtype=np.float64)
    orde = np.argsort(alles, kind="mergesort")
    gesorteerd = alles[orde]
    i = 0
    while i < len(gesorteerd):
        j = i
        while j + 1 < len(gesorteerd) and gesorteerd[j + 1] == gesorteerd[i]:
            j += 1
        r[orde[i:j + 1]] = (i + j) / 2.0 + 1.0
        i = j + 1
    nh = len(hoog)
    return (r[:nh].sum() - nh * (nh + 1) / 2.0) / (nh * len(laag))


def interval(hoog, laag, taxa_h, taxa_l, eenheid: str, n: int = 2000, zaad: int = 0):
    """Trek met teruglegging op de gevraagde eenheid en geef de 2,5e en 97,5e."""
    rng = np.random.default_rng(zaad)
    uit = []
    if eenheid == "rijen":
        for _ in range(n):
            a = hoog[rng.integers(0, len(hoog), len(hoog))]
            b = laag[rng.integers(0, len(laag), len(laag))]
            uit.append(auc(a, b))
    else:
        hs = {t: hoog[taxa_h == t] for t in np.unique(taxa_h)}
        ls = {t: laag[taxa_l == t] for t in np.unique(taxa_l)}
        hk, lk = np.array(list(hs)), np.array(list(ls))
        for _ in range(n):
            a = np.concatenate([hs[t] for t in hk[rng.integers(0, len(hk), len(hk))]])
            b = np.concatenate([ls[t] for t in lk[rng.integers(0, len(lk), len(lk))]])
            uit.append(auc(a, b))
    q = np.quantile(uit, [0.025, 0.975])
    return float(q[0]), float(q[1])


def groep_van(kop: str) -> str:
    for deel in kop.split():
        if deel.startswith("group="):
            return deel[6:]
    return "?"


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--crop", type=int, default=0,
                   help="elke reeks tot dit aantal basen (0 = laten staan)")
    p.add_argument("--nearest", action="store_true",
                   help="de naaste ENKELE reeks in plaats van de groepskern")
    a = p.parse_args()

    knip = (lambda s: s[:a.crop]) if a.crop else (lambda s: s)
    ref = lees("reference.fasta")
    groepen = [l.strip() for l in open("reference.groups")]
    if len(groepen) != len(ref):
        print("reference.groups past niet op reference.fasta", file=sys.stderr)
        return 1

    R = np.stack([profiel(knip(s)) for _, s in ref])
    if a.nearest:
        banken = R
    else:
        namen = sorted(set(groepen))
        g = np.array(groepen)
        banken = np.stack([R[g == n].mean(0) for n in namen])
        banken /= np.linalg.norm(banken, axis=1, keepdims=True)

    def scoor(pad):
        rijen = lees(pad)
        Q = np.stack([profiel(knip(s)) for _, s in rijen])
        return (Q @ banken.T).max(1), np.array([groep_van(k) for k, _ in rijen])

    sk, tk = scoor("known.fasta")
    sn, tn = scoor("novel.fasta")
    a_auc = auc(sk, sn)

    hoe = "naaste enkele reeks" if a.nearest else "een kern per groep"
    lengte = f"{a.crop} basen" if a.crop else "hele reeks"
    print(f"\n  dark-taxa-novelty v1.0   {hoe}, {lengte}")
    print(f"  {len(sk)} bekend tegen {len(sn)} nieuw, {len(set(groepen))} referentiegroepen")
    print(f"\n  AUC                      {a_auc:.4f}")
    for eenheid in ("rijen", "taxa"):
        lo, hi = interval(sk, sn, tk, tn, eenheid)
        print(f"  95% over {eenheid:<16}[{lo:.4f}, {hi:.4f}]")
    print(f"  gokken                   0.5000\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
