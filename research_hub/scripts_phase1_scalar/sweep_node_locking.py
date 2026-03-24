"""
Node-Locking Charge Hypothesis — Automated Sweep

Tests whether force direction correlates with node parity:
- ODD node separation → predict attraction (opposite charge)
- EVEN node separation → predict repulsion (same charge)

Sweeps WC separation from 1 to 30 nodes (0.5λ to 15λ),
reports MATCH/MISMATCH for each, and prints summary.
"""

import numpy as np
import wave_engine_1D_v3 as v3


def run_sweep(max_nodes=30):
    """Sweep node separations and test node-locking prediction."""
    results = []

    for sep_nodes in range(1, max_nodes + 1):
        # Position WCs symmetrically at node positions
        v3.wave_centers[0].x_am = v3._snap_to_node(-sep_nodes * v3.lam_am / 4)
        v3.wave_centers[1].x_am = v3._snap_to_node(+sep_nodes * v3.lam_am / 4)

        # Compute combined RMS and force
        rms = v3.compute_combined_rms(v3.x_am)
        force = v3.compute_force_field(rms) * v3.FORCE_TO_N

        # Force at each WC position
        idx_left = np.argmin(np.abs(v3.x_am - v3.wave_centers[0].x_am))
        idx_right = np.argmin(np.abs(v3.x_am - v3.wave_centers[1].x_am))
        f_left = force[idx_left]
        f_right = force[idx_right]

        # Determine actual force direction
        is_attract = f_left > 0 and f_right < 0
        is_repel = f_left < 0 and f_right > 0

        # Node-locking prediction
        parity = "ODD" if sep_nodes % 2 == 1 else "EVEN"
        expected = "attract" if sep_nodes % 2 == 1 else "repel"
        actual = "attract" if is_attract else "repel" if is_repel else "unclear"
        match = actual == expected

        # Node indices
        n1 = v3._get_node_index(v3.wave_centers[0].x_am)
        n2 = v3._get_node_index(v3.wave_centers[1].x_am)

        results.append({
            "sep_nodes": sep_nodes,
            "sep_lam": sep_nodes / 2,
            "parity": parity,
            "expected": expected,
            "actual": actual,
            "match": match,
            "f_left": f_left,
            "f_right": f_right,
            "n1": n1,
            "n2": n2,
        })

    return results


def print_results(results):
    """Print formatted results table."""
    print()
    print("=" * 90)
    print("NODE-LOCKING CHARGE HYPOTHESIS — SWEEP RESULTS")
    print("=" * 90)
    print(f"{'Sep':>4} {'λ':>6} {'Parity':>6} {'Expect':>8} {'Actual':>8} "
          f"{'Result':>10} {'F_left':>12} {'F_right':>12} {'Nodes':>10}")
    print("-" * 90)

    n_match = 0
    n_total = len(results)

    for r in results:
        tag = "✓ MATCH" if r["match"] else "✗ MISS"
        color_reset = ""
        print(f"{r['sep_nodes']:>4} {r['sep_lam']:>5.1f}λ {r['parity']:>6} "
              f"{r['expected']:>8} {r['actual']:>8} {tag:>10} "
              f"{r['f_left']:>12.3e} {r['f_right']:>12.3e} "
              f"  n{r['n1']}→n{r['n2']}")
        if r["match"]:
            n_match += 1

    print("-" * 90)
    print(f"TOTAL: {n_match}/{n_total} match ({100*n_match/n_total:.0f}%)")
    print()

    # Breakdown by parity
    odd_results = [r for r in results if r["parity"] == "ODD"]
    even_results = [r for r in results if r["parity"] == "EVEN"]
    odd_match = sum(1 for r in odd_results if r["match"])
    even_match = sum(1 for r in even_results if r["match"])

    print(f"ODD separations (expect attract):  {odd_match}/{len(odd_results)} match")
    print(f"EVEN separations (expect repel):   {even_match}/{len(even_results)} match")
    print()

    # Check if there's a consistent pattern even if not matching prediction
    odd_attract = sum(1 for r in odd_results if r["actual"] == "attract")
    odd_repel = sum(1 for r in odd_results if r["actual"] == "repel")
    even_attract = sum(1 for r in even_results if r["actual"] == "attract")
    even_repel = sum(1 for r in even_results if r["actual"] == "repel")

    print("Force direction distribution:")
    print(f"  ODD  seps: {odd_attract} attract, {odd_repel} repel, "
          f"{len(odd_results)-odd_attract-odd_repel} unclear")
    print(f"  EVEN seps: {even_attract} attract, {even_repel} repel, "
          f"{len(even_results)-even_attract-even_repel} unclear")
    print()

    # Is the pattern just alternating regardless of parity?
    actual_signs = []
    for r in results:
        if r["actual"] == "attract":
            actual_signs.append(+1)
        elif r["actual"] == "repel":
            actual_signs.append(-1)
        else:
            actual_signs.append(0)

    alternating = 0
    for i in range(1, len(actual_signs)):
        if actual_signs[i] != 0 and actual_signs[i-1] != 0:
            if actual_signs[i] != actual_signs[i-1]:
                alternating += 1

    print(f"Force alternates between consecutive separations: "
          f"{alternating}/{n_total-1} transitions")

    if alternating > 0.8 * (n_total - 1):
        print("→ Force DOES alternate every node step (λ/2)")
        if n_match > 0.8 * n_total:
            print("→ AND matches node-locking prediction — HYPOTHESIS SUPPORTED")
        else:
            print("→ BUT phase doesn't match even/odd prediction — "
                  "alternation exists but parity assignment may be shifted")
    else:
        print("→ Force does NOT consistently alternate — more complex pattern")

    print("=" * 90)


if __name__ == "__main__":
    print("Running node-locking sweep...")
    results = run_sweep(30)
    print_results(results)
