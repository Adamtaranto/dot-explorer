"""Measure k-mer index memory for synthetic assembly pairs.

Answers "where does k-mer mode memory go, and how large an input fits a
given heap?".  For each requested size the script forks a clean worker
subprocess that builds a ``CrossIndex`` over a homologous query/target pair
(the benchmarks' generator), computes stranded matches, and reports:

- peak RSS of the whole build (``ru_maxrss``, i.e. what the OS saw),
- final RSS after the build,
- the index's own accounting from ``CrossIndex.approx_bytes()``
  (``seq_bytes`` / ``kmer_index`` / ``pair_cache``),
- the number of match records materialised.

Run from the repo root with the compiled extension importable::

    python scripts/mem_profile_index.py                 # 10, 50, 100 Mbp
    python scripts/mem_profile_index.py --mb 10 50      # subset
    python scripts/mem_profile_index.py --min-block-len 50

Bytes-per-bp from the table extrapolate to the browser app's wasm-heap
budget (see ``_KMER_HARD_LIMIT`` in app/app.py).
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import resource
import subprocess
import sys

SIZES_MB = (10, 50, 100)
K = 15  # the app's default k-mer length


def _rss_bytes(ru_maxrss: int) -> int:
    """Normalise ``ru_maxrss`` to bytes (macOS reports bytes, Linux KiB)."""
    return ru_maxrss if sys.platform == 'darwin' else ru_maxrss * 1024


def _current_rss() -> int:
    """Return the process's current resident set size in bytes."""
    ru = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return _rss_bytes(ru)


def worker(length: int, k: int, min_block_len: int) -> dict:
    """Build one index pair in this process and return its measurements.

    Parameters
    ----------
    length : int
        Per-sequence length in bases (query and target each).
    k : int
        K-mer length.
    min_block_len : int
        Native match-block length filter passed to ``compute_matches``.

    Returns
    -------
    dict
        Measurements; see module docstring.
    """
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'python/benchmarks'))
    from _synth import homologous_pair

    from rusty_dot import CrossIndex

    query, target = homologous_pair(length, seed=42, divergence=0.02)
    baseline = _current_rss()

    index = CrossIndex(k)
    index.add_sequence('q1', query, group='query')
    index.add_sequence('t1', target, group='target')
    # Release the Python strings: Rust holds its own copy (this mirrors the
    # app's provider-driven build, which never retains parsed strings).
    del query, target

    index.compute_matches(
        query_group='query',
        target_group='target',
        merge=True,
        min_block_len=min_block_len,
    )
    n_records = len(index.get_records_for_pair('query', 'target'))

    peak = _rss_bytes(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return {
        'length': length,
        'baseline_rss': baseline,
        'peak_rss': peak,
        'approx': index.approx_bytes(),
        'n_records': n_records,
    }


def main() -> None:
    """Orchestrate one clean worker subprocess per size and print a table."""
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--mb', type=int, nargs='+', default=list(SIZES_MB))
    ap.add_argument('--k', type=int, default=K)
    ap.add_argument('--min-block-len', type=int, default=0)
    ap.add_argument('--worker', type=int, help=argparse.SUPPRESS)
    args = ap.parse_args()

    if args.worker is not None:
        print(json.dumps(worker(args.worker, args.k, args.min_block_len)))
        return

    header = (
        f'{"size":>8} {"peak RSS":>10} {"peak B/bp":>9} {"index":>10} '
        f'{"seq":>8} {"kmer":>8} {"pairs":>8} {"records":>9}'
    )
    print(f'k={args.k}, min_block_len={args.min_block_len} (per-side sizes)')
    print(header)
    for mb in args.mb:
        length = mb * 1_000_000
        out = subprocess.run(
            [
                sys.executable,
                __file__,
                '--worker',
                str(length),
                '--k',
                str(args.k),
                '--min-block-len',
                str(args.min_block_len),
            ],
            capture_output=True,
            text=True,
            check=True,
        )
        r = json.loads(out.stdout.strip().splitlines()[-1])
        total_bp = 2 * length
        approx = r['approx']
        gb = 1024**3
        mby = 1024**2
        print(
            f'{mb:>5} Mb {r["peak_rss"] / gb:>8.2f}G '
            f'{r["peak_rss"] / total_bp:>9.1f} '
            f'{approx["total"] / gb:>9.2f}G '
            f'{approx["seq_bytes"] / mby:>7.0f}M '
            f'{approx["kmer_index"] / mby:>7.0f}M '
            f'{approx["pair_cache"] / mby:>7.0f}M '
            f'{r["n_records"]:>9,}'
        )


if __name__ == '__main__':
    main()
