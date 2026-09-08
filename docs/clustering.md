# Similarity & Clustering

rusty-dot can compare every sequence in a set against every other with
[sourmash](https://sourmash.readthedocs.io) FracMinHash sketches, cluster
them hierarchically, draw the resulting tree beside the dot-plot matrix or
a similarity heatmap, and assign contigs to clusters at a cutoff you
choose. Install the optional dependencies with:

```bash
pip install "rusty-dot[cluster]"    # sourmash + scipy
```

This page explains how to choose between the similarity metrics and
clustering modes; the [clustering tutorial](tutorials/clustering.md) walks
through the full workflow in code.

## How sketching works

Every metric starts from a **FracMinHash sketch** of each sequence
(`compute_sketches`): the sequence's canonical k-mers are hashed, and
roughly one in every `scaled` hashes is kept. Three parameters matter:

- **`ksize`** (default 21) — the k-mer length. 21 is the sourmash default
  for DNA and discriminates well between genomes; use 31 for
  strain-level comparisons of close genomes and smaller k only when
  comparing highly diverged sequences.
- **`scaled`** (default 1000) — the sampling rate. One hash per ~1000 bp
  is fine for genome-scale comparisons; for short sequences, or when you
  need precise ANI estimates, use a **lower** value (e.g. 100, or even 10
  for sequences of a few kb) so each sketch keeps enough hashes to
  compare. As a rule of thumb you want at least a few hundred hashes per
  sketch (`len(seq) / scaled`).
- **`track_abundance`** (default on, sourmash's `-p abund`) — record how
  often each kept k-mer occurs. Required for angular similarity, ignored
  by everything else, and cheap to keep on.

## Choosing a similarity metric

| Metric | What it measures | Needs abundance | Symmetric |
|---|---|---|---|
| `jaccard` | Fraction of k-mers shared between the two sequences (like `sourmash compare`) | no | yes |
| `angular` | Cosine-style similarity weighted by k-mer abundance (like `sourmash compare` with abundance) | **yes** | yes |
| `ani` | Estimated average nucleotide identity, from max containment, with 95% confidence bounds | no | yes |
| `containment` | Fraction of sequence *i*'s k-mers found in *j* (like `sourmash gather`'s containment) | no | **no** |
| `max_containment` / `avg_containment` | Symmetrised containment (max/mean of the two directions) | no | yes |

Guidance:

- **Jaccard** is the default general-purpose choice: it treats the two
  sequences symmetrically and penalises both unshared content and size
  differences. Two sequences where one is a perfect subset of the other
  score *low* Jaccard — if that is not what you want, use containment.
- **Angular similarity** behaves like Jaccard but weights repeated
  k-mers by how often they occur, so it distinguishes sequences that
  share k-mer *sets* but differ in repeat copy number (e.g. satellite
  arrays, transposon load).
- **ANI** approximates percent identity, the quantity most people mean
  by "80% identical". It is derived from containment under a Poisson
  model, and rusty-dot always reports the 95% confidence interval
  alongside — **use a lower `scaled` to tighten it**. Pairs that share
  no hashes get ANI 0 (the estimator is undefined there).
- **Containment** is the k-mer analogue of *percentage coverage*: how
  much of sequence A is present in B, regardless of how much extra
  content B carries. It is asymmetric — a plasmid is fully contained in
  the genome that carries it, not vice versa — which is exactly what the
  dual identity + coverage clustering mode (below) exploits.

## Clustering

### Hierarchical tree + similarity cutoff

`linkage_from_similarity` converts a similarity matrix to distances
(`1 − similarity`) and runs scipy hierarchical clustering (UPGMA /
`average` linkage by default). `Tree.from_linkage` turns the linkage
into a tree you can pass to `DotPlotter.plot(tree=...)` or
`plot_similarity_heatmap(tree=...)`, and `assign_clusters(sim, cutoff)`
cuts it: contigs whose linkage-aggregated similarity is at least
`cutoff` share a cluster. Pass `tree_cutoff=1 - cutoff` to draw the cut
as a dashed line through the dendrogram.

Because the tree and the assignment come from the same linkage, cluster
members are always contiguous along the tree order, so cluster borders
form clean blocks in the matrix.

### Dual identity + coverage thresholds

`assign_clusters_dual(identity, coverage, ...)` links two sequences only
when **both** conditions hold — for example ANI ≥ 80% *and* containment
≥ 80% — and clusters are the connected components of the link graph.
This matches the common "80/80" rule for grouping related elements:

```python
ani = pairwise_similarity(sketches, metric='ani')
cov = pairwise_similarity(sketches, metric='containment')
clusters = assign_clusters_dual(
    ani, cov, identity_cutoff=0.80, coverage_cutoff=0.80, reciprocal=True
)
```

`reciprocal=True` (default) requires the coverage threshold in *both*
directions — a short fragment nested in a long sequence does **not**
join its cluster, because the long sequence is poorly covered by the
fragment. With `reciprocal=False` one passing direction suffices, which
pulls nested fragments into their parents' clusters.

**Caveat — containment is exact-k-mer coverage.** Uniformly scattered
SNPs destroy k-mers steeply: at 98% identity, only about
`0.98²¹ ≈ 65%` of 21-mers survive intact, so a fully-overlapping but
SNP-diverged pair shows ~65% containment, not ~100%. Expected
containment is roughly `coverage × ANI^k`. The 80/80 rule with sketch
containment therefore groups sequences that share long (near-)exact
stretches — structural variants, redundant haplotigs, nested elements.
For SNP-diverged families either lower the coverage cutoff toward
`0.8 × ANI^k`, or cluster on a plain similarity cutoff instead; for a
strict alignment-based 80/80 (CheckV-style ANI + aligned fraction),
compute ANI outside rusty-dot from real alignments.

Note that connected components chain: A links B, B links C ⇒ {A, B, C}
share a cluster even if A and C do not pass the thresholds directly.
Dual-mode clusters are also not guaranteed to be contiguous along a
tree's leaf order; non-contiguous clusters are outlined per block, with a
warning.

## Putting it on a plot

- `DotPlotter.plot(tree=...)` draws the dendrogram left of the matrix
  rows and fixes the row (and, for self-comparisons, column) order to
  the tree — `contig_order` and `auto_reverse` are rejected while a tree
  is active, and the contig name labels move onto the tree tips.
  `cluster_borders=` takes a `ClusterResult` and outlines each cluster's
  block of panels.
- `plot_similarity_heatmap(sim, tree=..., clusters=...)` renders the
  matrix itself as a heatmap with the tree at the left, cluster
  outlines, a selectable colormap, and a scale bar at the right.
- User-supplied trees (newick or IQ-TREE `.treefile`, via `Tree.read`)
  work the same way; tip labels must match the sequence names exactly
  (mismatches raise an error listing both directions).
- `SimilarityMatrix.to_csv` and `ClusterResult.to_csv` export the raw
  matrix and the contig → cluster table.

See the [API reference](api/similarity.md) for full signatures and the
[tutorial](tutorials/clustering.md) for a worked example.
