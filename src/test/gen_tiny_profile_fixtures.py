#!/usr/bin/env python3
"""Deterministic generator for the tiny profile-filter test fixtures.

Regenerates the three checked-in fixtures used by
TestTinyProfileCpuFilter.sh / TestTinyProfileGpuFilter.sh:

  tiny_mix.sto      10 x 60 Stockholm, residues drawn from the BLOSUM62
                    background (seed 20260919).
  tiny_polyW.sto    10 x 60 Stockholm, ~90% W + 10% background (seed 77).
  tiny_target.fasta 18 random targets, lengths cycling 60/80/100 (seed 7),
                    plus spike_exact = tiny_mix row 0 and spike_mut70 =
                    70% identity to spike_exact (18 substitutions, seed 99).

The BLOSUM62 background is read from data/blosum62.out (line 2), so there are
no duplicated constants. Output is byte-for-byte reproducible on a fixed
CPython, because it uses only random.Random(seed).random(); it is NOT stable
across Python versions, so record the version below whenever the fixtures are
regenerated.

Generated with: Python 3.12.3

Usage:
  python3 gen_tiny_profile_fixtures.py [outdir]          # default: script dir
  python3 gen_tiny_profile_fixtures.py --check [outdir]  # byte-compare only
"""
import os
import random
import sys

AAS = "ACDEFGHIKLMNPQRSTVWY"
MIX_SEED = 20260919
POLYW_SEED = 77
TARGET_SEED = 7
SPIKE_SEED = 99
MIX_NSEQ = 10
MIX_LEN = 60
POLYW_NSEQ = 10
POLYW_LEN = 60
POLYW_W_FRAC = 0.90
TARGET_LENS = [60, 80, 100]
SPIKE_IDENTITY = 0.70

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))


def load_background():
    path = os.path.join(ROOT, "data", "blosum62.out")
    with open(path) as fh:
        for line in fh:
            if line.startswith("# Background"):
                nums = [float(x) for x in line.split(":")[1].split()]
                assert len(nums) == len(AAS) + 1, "unexpected background size"
                probs = nums[:len(AAS)]  # drop the trailing X
                total = sum(probs)
                assert total > 0.0, "empty background"
                return probs, total
    raise RuntimeError("no '# Background' line in %s" % path)


def _pick(rng, probs, total):
    u = rng.random() * total
    acc = 0.0
    for i, p in enumerate(probs):
        acc += p
        if u < acc:
            return AAS[i]
    return AAS[-1]


def _random_seq(rng, probs, total, length):
    return "".join(_pick(rng, probs, total) for _ in range(length))


def _mutate(rng, probs, total, seq, nMut):
    chars = list(seq)
    positions = list(range(len(chars)))
    rng.shuffle(positions)
    for pos in positions[:nMut]:
        alt = chars[pos]
        while alt == chars[pos]:
            alt = _pick(rng, probs, total)
        chars[pos] = alt
    return "".join(chars)


def build():
    probs, total = load_background()

    mix_rng = random.Random(MIX_SEED)
    mix = [_random_seq(mix_rng, probs, total, MIX_LEN) for _ in range(MIX_NSEQ)]

    pw_rng = random.Random(POLYW_SEED)
    polyw = []
    for _ in range(POLYW_NSEQ):
        row = []
        for _ in range(POLYW_LEN):
            row.append("W" if pw_rng.random() < POLYW_W_FRAC
                       else _pick(pw_rng, probs, total))
        polyw.append("".join(row))

    tgt_rng = random.Random(TARGET_SEED)
    targets = []
    for i in range(18):
        length = TARGET_LENS[i % len(TARGET_LENS)]
        targets.append(("t%d" % i, _random_seq(tgt_rng, probs, total, length)))
    spike_exact = mix[0]
    spike_rng = random.Random(SPIKE_SEED)
    nMut = int(round((1.0 - SPIKE_IDENTITY) * len(spike_exact)))
    targets.append(("spike_exact", spike_exact))
    targets.append(("spike_mut70", _mutate(spike_rng, probs, total,
                                           spike_exact, nMut)))

    mix_sto = ["# STOCKHOLM 1.0", "#=GF ID tiny_mix"]
    mix_sto += ["mix%d %s" % (i, s) for i, s in enumerate(mix)]
    mix_sto.append("//")

    pw_sto = ["# STOCKHOLM 1.0", "#=GF ID tiny_polyW"]
    pw_sto += ["pw%d %s" % (i, s) for i, s in enumerate(polyw)]
    pw_sto.append("//")

    fasta = []
    for name, seq in targets:
        fasta.append(">%s" % name)
        fasta.append(seq)

    return {
        "tiny_mix.sto": "\n".join(mix_sto) + "\n",
        "tiny_polyW.sto": "\n".join(pw_sto) + "\n",
        "tiny_target.fasta": "\n".join(fasta) + "\n",
    }


def main(argv):
    check = False
    args = argv[1:]
    if args and args[0] == "--check":
        check = True
        args = args[1:]
    outdir = args[0] if args else HERE
    files = build()
    status = 0
    for name in sorted(files):
        content = files[name]
        path = os.path.join(outdir, name)
        if check:
            try:
                with open(path) as fh:
                    existing = fh.read()
            except OSError as exc:
                print("MISSING %s (%s)" % (path, exc))
                status = 1
                continue
            if existing == content:
                print("ok %s" % path)
            else:
                print("DIFF %s" % path)
                status = 1
        else:
            with open(path, "w") as fh:
                fh.write(content)
            print("wrote %s" % path)
    if not check:
        print("generated with Python %s" % sys.version.split()[0])
    return status


if __name__ == "__main__":
    sys.exit(main(sys.argv))
