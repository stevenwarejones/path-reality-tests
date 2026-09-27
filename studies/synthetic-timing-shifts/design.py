#!/usr/bin/env python3
"""Synthetic sparse-record injections; no NIST events are read or modified."""
import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.special import ndtr
from scipy.stats import beta, binom

HERE = Path(__file__).resolve().parent
TESTS = ("click_count", "pooled_timing", "trained_timing", "trained_histogram")
# All real arrival times have an outcome, including both overflow tails.
EDGES = np.r_[-np.inf, np.arange(-8, 9, dtype=float), np.inf]
TIMING_SIGN = np.where(EDGES[1:] <= 0, -1, 1).astype(np.int64)


@dataclass(frozen=True)
class Config:
    trials: int = 53_554_798  # half the archive size; illustrative fixed receiver setting
    blocks: int = 32
    training_blocks: int = 8
    click_probability: float = 0.0002
    timing_sigma_bins: float = 2.0
    tag_bin_ps: float = 78.125
    family_alpha: float = 0.01

    def __post_init__(self):
        if not isinstance(self.trials, int) or self.trials < self.blocks:
            raise ValueError("trials must be an integer at least blocks")
        if (not isinstance(self.blocks, int) or not isinstance(self.training_blocks, int)
                or self.blocks < 4 or self.blocks % 2 or self.training_blocks % 2
                or not 0 < self.training_blocks < self.blocks):
            raise ValueError("even block counts and nonempty training/test segments required")
        if not 0 < self.click_probability < 0.5:
            raise ValueError("click_probability must be between 0 and 0.5")
        if not np.isfinite(self.timing_sigma_bins) or self.timing_sigma_bins <= 0:
            raise ValueError("positive finite timing sigma required")
        if not np.isfinite(self.tag_bin_ps) or self.tag_bin_ps <= 0:
            raise ValueError("positive finite time-bin scale required")
        if not 0 < self.family_alpha < 1:
            raise ValueError("family_alpha must lie in (0, 1)")


@dataclass(frozen=True)
class Scenario:
    name: str
    shift_bins: float = 0.0  # difference between the two assignment-arm means
    rate_gap_fraction: float = 0.0  # (p1-p0)/baseline, before reversal
    reverse: bool = False
    drift: bool = False
    transfer_failure: bool = False

    def __post_init__(self):
        if not np.isfinite(self.shift_bins):
            raise ValueError("finite shift required")
        if not np.isfinite(self.rate_gap_fraction) or abs(self.rate_gap_fraction) > 1:
            raise ValueError("rate gap fraction must be in [-1, 1]")


def timing_probabilities(mean, sigma):
    if not np.isfinite(mean) or not np.isfinite(sigma) or sigma <= 0:
        raise ValueError("finite mean and positive finite sigma required")
    p = np.diff(ndtr((EDGES - mean) / sigma))
    return p / p.sum()


def cell_law(config, scenario, block, setting):
    """Return click probability and unconditional-on-window arrival distribution."""
    state = block % 2  # observable pretrial regime, not learned from current outcomes
    direction = (1 if state == 0 else -1) if scenario.reverse else 1
    # Deliberate counterexample to transport: same phase in training, opposite
    # phases across the two regimes at test. An unsigned two-sided score cannot
    # fix this loss of association with the learned feature.
    if scenario.transfer_failure:
        direction = 1 if block < config.training_blocks or state == 0 else -1
    t = block / (config.blocks - 1)
    baseline = config.click_probability * (0.5 + t if scenario.drift else 1)
    common_clock = 3 * (2 * t - 1) if scenario.drift else 0
    arm = 2 * setting - 1
    rate = baseline * (1 + arm * direction * scenario.rate_gap_fraction / 2)
    mean = common_clock + arm * direction * scenario.shift_bins / 2
    return rate, timing_probabilities(mean, config.timing_sigma_bins)


def simulate(config, scenario, rng):
    """Exact count compression of independent trials, not a Poisson approximation.

    Each trial: fair setting, Bernoulli click, Gaussian arrival-time bin if clicked.
    No individual event arrays or 100-million-row allocation are necessary.
    """
    sizes = np.full(config.blocks, config.trials // config.blocks, dtype=np.int64)
    sizes[:config.trials % config.blocks] += 1
    trials = np.empty((config.blocks, 2), dtype=np.int64)
    clicks = np.zeros((config.blocks, 2, len(TIMING_SIGN)), dtype=np.int64)
    for block, n in enumerate(sizes):
        trials[block, 0] = rng.binomial(n, 0.5)
        trials[block, 1] = n - trials[block, 0]
        for setting in (0, 1):
            rate, time_law = cell_law(config, scenario, block, setting)
            k = rng.binomial(trials[block, setting], rate)
            clicks[block, setting] = rng.multinomial(k, time_law)
    return trials, clicks


def train_features(trials, clicks):
    """Only the initial segment enters this API; returns frozen integer weights."""
    weights = np.zeros((2, clicks.shape[-1]), dtype=np.int64)
    timing = np.zeros_like(weights)
    for state in (0, 1):
        counts = clicks[state::2].sum(axis=0)
        n = trials[state::2].sum(axis=0)
        if np.any(n == 0):
            continue
        difference = counts[1] / n[1] - counts[0] / n[0]
        weights[state] = np.sign(difference).astype(np.int64)
        timing[state] = int(np.sign(difference @ TIMING_SIGN)) * TIMING_SIGN
    return {"trained_timing": timing, "trained_histogram": weights}


def exact_score_pvalue(positive, negative):
    """Two-sided fair-sign test, conditional on the number of nonzero terms."""
    if any(not isinstance(v, (int, np.integer)) or v < 0 for v in (positive, negative)):
        raise ValueError("nonnegative integer score counts required")
    n = int(positive + negative)
    return 1.0 if n == 0 else min(1.0, float(2 * binom.cdf(min(positive, negative), n, 0.5)))


def score(clicks, weights, first_block):
    positive = negative = 0
    for local_block, counts in enumerate(clicks):
        w = weights[(first_block + local_block) % 2]
        positive += int(counts[1, w > 0].sum() + counts[0, w < 0].sum())
        negative += int(counts[0, w > 0].sum() + counts[1, w < 0].sum())
    return exact_score_pvalue(positive, negative)


def evaluate(config, trials, clicks):
    split = config.training_blocks
    features = train_features(trials[:split], clicks[:split])
    features["click_count"] = np.ones((2, len(TIMING_SIGN)), dtype=np.int64)
    features["pooled_timing"] = np.tile(TIMING_SIGN, (2, 1))
    # Every detector gets the same held-out segment. No-clicks have score zero.
    return {name: score(clicks[split:], features[name], split) for name in TESTS}


def binomial_interval(successes, repeats):
    """Pointwise 95% Clopper-Pearson Monte Carlo interval (not a physics bound)."""
    return [0.0 if successes == 0 else float(beta.ppf(.025, successes, repeats-successes+1)),
            1.0 if successes == repeats else float(beta.ppf(.975, successes+1, repeats-successes))]


def scenarios():
    cases = [Scenario("null"), Scenario("null_with_drift", drift=True)]
    for reverse in (False, True):
        for shift in (0.125, 0.25, 0.5, 1.0):
            cases.append(Scenario(f"{'reversing' if reverse else 'fixed'}_shift_{shift:g}",
                                  shift_bins=shift, reverse=reverse))
    cases.extend([
        Scenario("fixed_rate_gap", rate_gap_fraction=0.5),
        Scenario("reversing_rate_gap", rate_gap_fraction=0.5, reverse=True),
        Scenario("regime_transfer_failure", shift_bins=1, transfer_failure=True),
    ])
    return cases


def run(repeats=400, seed=20260926):
    if not isinstance(repeats, int) or repeats <= 0:
        raise ValueError("positive integer repeat count required")
    rows = []
    for n in (10_000_000, 53_554_798):
        config = Config(trials=n)
        for case in scenarios():
            # Stable per-scenario seed, unaffected by adding/reordering other cases.
            label = f"{seed}:{n}:{case.name}"
            entropy = int.from_bytes(hashlib.sha256(label.encode()).digest()[:16], "little")
            rng = np.random.Generator(np.random.PCG64(entropy))
            rejected = dict.fromkeys(TESTS, 0)
            family = train_clicks = test_clicks = 0
            for _ in range(repeats):
                trials, clicks = simulate(config, case, rng)
                pvalues = evaluate(config, trials, clicks)
                decisions = {name: p <= config.family_alpha / len(TESTS)
                             for name, p in pvalues.items()}
                for name, decision in decisions.items():
                    rejected[name] += int(decision)
                family += int(any(decisions.values()))
                train_clicks += int(clicks[:config.training_blocks].sum())
                test_clicks += int(clicks[config.training_blocks:].sum())
            rows.append({"config": asdict(config), "scenario": asdict(case),
                         "repeats": repeats, "rejections": rejected,
                         "rejection_intervals_95": {name: binomial_interval(k, repeats)
                                                    for name, k in rejected.items()},
                         "family_rejections": family,
                         "family_interval_95": binomial_interval(family, repeats),
                         "mean_training_clicks": train_clicks / repeats,
                         "mean_test_clicks": test_clicks / repeats})
    return {"schema_version": 1, "data_kind": "synthetic_only", "seed": seed,
            "rng": "NumPy PCG64; SHA256-derived per-scenario seeds",
            "experimental_claim_enabled": False, "causal_claim_enabled": False,
            "spacelike_claim_enabled": False, "rows": rows}


def render_report(result):
    lines = ["# Artificial timing shifts: synthetic sensitivity", "",
             "Generated by `design.py`. These are Monte Carlo detection frequencies,",
             "not measured NIST effects or apparatus-certified limits. See [design and",
             "assumptions](README.md) before interpreting the picosecond scale.", "",
             "Each detector uses the final 24 of 32 blocks; the first 8 train the two",
             "learned scores. Four two-sided tests share family alpha 0.01 by Bonferroni.",
             "Timing shifts are the difference between assignment-arm means; reversals",
             "alternate across an observable pretrial regime. No-clicks remain in the",
             "population. All timing bins, including overflow tails, are retained.", "",
             "The assumed click rate is 200 per million and Gaussian width is 2 tag bins",
             "(156.25 ps). Neither width nor regime persistence is calibrated to NIST.", ""]
    for n in sorted({row["config"]["trials"] for row in result["rows"]}):
        rows = [r for r in result["rows"] if r["config"]["trials"] == n]
        null = rows[0]
        lines.extend([f"## {n:,} synthetic trials", "",
                      f"{null['repeats']} independent simulations per row. Mean clicks under the stationary null: "
                      f"{null['mean_training_clicks']:.1f} training, {null['mean_test_clicks']:.1f} test.", "",
                      "| Injection | Shift (ps) | Count | Pooled timing | Trained timing | Trained histogram | Any test (95% MC interval) |",
                      "|---|---:|---:|---:|---:|---:|---|",
                      ])
        for row in rows:
            case = row["scenario"]
            rates = [100 * row["rejections"][name] / row["repeats"] for name in TESTS]
            family = 100 * row["family_rejections"] / row["repeats"]
            lo, hi = [100 * p for p in row["family_interval_95"]]
            shift = case["shift_bins"] * row["config"]["tag_bin_ps"]
            lines.append(f"| {case['name']} | {shift:.2f} | " +
                         " | ".join(f"{p:.1f}%" for p in rates) +
                         f" | {family:.1f}% ({lo:.1f}–{hi:.1f}%) |")
        lines.append("")
    lines.extend(["## Interpretation limits", "",
                  "Null rows measure false alarms. Other rows measure power for the exact",
                  "injected model. Rate-gap injections use p1-p0 = 0.5 times the baseline",
                  "(100 clicks per million), with opposite signs in the reversing case.",
                  "The transfer-failure row intentionally breaks the relationship between",
                  "the training feature and the test regimes; low power there is expected.", "",
                  "Intervals are pointwise 95% Monte Carlo intervals; all individual-test",
                  "intervals and integer rejection counts are in [results.json](results.json).",
                  "These are not simultaneous confidence bands across the scenario grid.",
                  "The family error bound applies separately to each synthetic experiment,",
                  "not to a scan of all shifts on one real dataset. A non-rejection does not",
                  "certify absence of dependence. No experimental record was tested here.", ""])
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--output-dir", type=Path, default=HERE)
    args = parser.parse_args()
    result = run()
    outputs = {"results.json": json.dumps(result, indent=2, sort_keys=True) + "\n",
               "report.md": render_report(result)}
    if args.check:
        for name, text in outputs.items():
            if (args.output_dir / name).read_text() != text:
                raise SystemExit(f"Stale synthetic artifact: {name}")
        print(f"Reproduced {len(result['rows'])} synthetic scenarios and report")
    else:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        for name, text in outputs.items():
            (args.output_dir / name).write_text(text)
        print(f"Wrote {len(result['rows'])} synthetic scenarios to {args.output_dir}")


if __name__ == "__main__":
    main()
