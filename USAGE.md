# Usage

Command-line usage for the five executables built by this project (see
`README.md` for build instructions). This reflects what each program's
argument parser actually accepts, not just the text of its `-h`/usage
message — a few programs have parsers that are stricter (or looser) than
their printed usage string, and that's noted inline.

---

## RunSimulation

```
RunSimulation -BCUYclhpt -n <ntax> -s <nsites> -r <height> -R <replicates> [-i <iterations>] [-P <rate>] [-w <width>]
```

Simulates trees and sequence alignments, then compares MCMC branch-length
sampling under the pruning algorithm vs. the LvD (decomposition) algorithm.

The **first argument** must be a single `-`-prefixed flag string (flags are
matched by character, order doesn't matter, unrecognised characters are
silently ignored):

| Flag | Meaning |
|------|---------|
| `B` | Simulate from beta-critical distribution |
| `C` | Simulate from random caterpillar |
| `U` | Simulate from uniform distribution |
| `Y` | Simulate from Yule distribution |
| `l` | Run/time the MCMC using the LvD (decomposition) algorithm |
| `p` | Run/time the MCMC using the pruning algorithm |
| `h` | Print a header row before results |
| `t` | Also print simulated trees (Newick) to stderr |

At least one of `B`, `C`, `U`, `Y` should be set, or no replicates are run.

Remaining arguments are `key value` pairs, in any order:

| Option | Value | Required | Default |
|--------|-------|----------|---------|
| `-n <ntax>` | integer, number of taxa | yes | — |
| `-s <nsites>` | integer, number of sites | yes | — |
| `-r <height>` | double, expected root-to-tip mutations | yes | — |
| `-R <replicates>` | integer, number of replicates | yes | — |
| `-i <iterations>` | integer, MCMC iterations | no | 1000 |
| `-P <rate>` | double, exponential prior rate on branch lengths | no | 20.0 |
| `-w <width>` | double > 0, half-width of uniform branch-length proposal | no | 0.01 |

Example:

```
RunSimulation -UYlph -n 100 -s 1000 -r 0.05 -R 20 -i 2000
```


---

## TreeHeightSimulation

```
TreeHeightSimulation -BCUYh -m <min_ntax> -n <max_ntax> -R <replicates>
```

Simulates trees only (no sequences/MCMC) and reports pruning-tree vs.
LvD-tree heights, sampling the number of taxa per replicate so that
`log(ntax)` is uniform between `min_ntax` and `max_ntax`.

First argument is the flag string:

| Flag | Meaning |
|------|---------|
| `B` | Simulate from beta-critical distribution |
| `C` | Simulate from random caterpillar |
| `U` | Simulate from uniform distribution |
| `Y` | Simulate from Yule distribution |
| `h` | Print a header row before results |

Remaining `key value` pairs:

| Option | Value | Required | Default |
|--------|-------|----------|---------|
| `-m <min_ntax>` | integer ≥ 1 | yes | 10 |
| `-n <max_ntax>` | integer, must be ≥ `min_ntax` | yes | 10000 |
| `-R <replicates>` | integer, number of replicates | yes | — |

Example:

```
TreeHeightSimulation -UYh -m 10 -n 5000 -R 50
```

---

## TreeHeights

```
TreeHeights [-h] <treefile> [<treefile> ...]
```

Reads one or more Newick tree files and reports pruning-tree vs. LvD-tree
heights for every tree in each file.

- `-h`, if present, **must be the first argument**; it enables a header row
  and is otherwise consumed (not treated as a filename).
- All remaining arguments are treated as tree filenames; any file that
  fails to open is skipped with a warning to stderr.
- Requires at least one argument (a `-h` alone with no filenames does not
  error, but processes zero files).

Example:

```
TreeHeights -h tree1.nwk tree2.nwk tree3.nwk
```

---

## LvDLikelihoods

```
LvDLikelihoods [-h] <seqfile> <treefile>
```

Reads a PHYLIP sequence file and a Newick tree file, and reports
likelihoods and timings under the pruning algorithm vs. the LvD algorithm
for every tree in the file.

- Accepts exactly 2 positional arguments (`<seqfile> <treefile>`), or 3
  arguments where the first is literally `-h` followed by the same two
  filenames.
- The header row is printed only when `-h` is given; it is omitted in the
  2-argument form.
- Any other argument count is rejected.

Examples:

```
LvDLikelihoods alignment.phy trees.nwk
LvDLikelihoods -h alignment.phy trees.nwk
```

---

## BranchLengthMCMC

```
BranchLengthMCMC -SPLh -p <prior> -u <proposal> -i <iterations> <seqfile> <treefile>
```

Runs Metropolis-Hastings MCMC on branch lengths of the first tree in
`<treefile>`, comparing up to three likelihood implementations: standard
pruning, decomposition-with-pruning-layout, and decomposition-with-LvD-
layout.

This parser is strict: it requires **exactly 10** arguments in the exact
positions below (not free-form `key value` pairs like `RunSimulation`).

| Position | Value | Notes |
|----------|-------|-------|
| 1 | flag string, e.g. `-SPLh` | must start with `-` |
| 2 | literal `-p` | |
| 3 | `<prior>` | double > 0, exponential prior rate λ (no default — always required) |
| 4 | literal `-u` | |
| 5 | `<proposal>` | double > 0, half-width of uniform proposal |
| 6 | literal `-i` | |
| 7 | `<iterations>` | integer > 0, number of MCMC iterations |
| 8 | `<seqfile>` | PHYLIP sequence file |
| 9 | `<treefile>` | Newick tree file |

Flags (characters in argument 1):

| Flag | Meaning |
|------|---------|
| `S` | Run MCMC using the standard likelihood |
| `P` | Run MCMC using the decomposition tree, pruning layout |
| `L` | Run MCMC using the decomposition tree, LvD layout |
| `h` | Print a header row |

At least one of `S`, `P`, `L` must be set.

Example:

```
BranchLengthMCMC -SPLh -p 20 -u 0.01 -i 5000 alignment.phy tree.nwk
```
