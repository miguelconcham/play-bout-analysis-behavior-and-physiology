# Behavior probability by target actor

This is an independent estimate/analyse pipeline. It does not overwrite or
change the original pipeline in `Figure 1/Behavior probability`.

## Data organization

Run `Estimate_actor_behavior_probability.m`.

Within every `Session_structure / Partners session`, the estimate finds every
unique real animal name. It creates one `actor_data` block per animal, including
animals that did not perform a selected target behavior.

For target actor A:

- `self_behavior`: behavior performed by A.
- `other_behavior`: behavior performed by any animal other than A. Multiple
  other animals are combined with a logical OR.

Both arrays have dimensions:

```text
target events performed by A × behavior types × time bins
```

The corresponding `*_behavior_wrapped` arrays contain a real-time pre-onset
flank, a fixed-bin representation from target onset through target offset, and
a real-time post-offset flank. They support pre-onset baseline comparisons
without treating the target duration as baseline.

The estimate is saved separately as:

```text
Data/Analysis results/behavior probability by target actor/
    actor_behavior_probability_around_offset.mat
```

## Analysis

Run `Analyze_actor_behavior_probability.m`.

It produces:

- Self and other heatmaps for every target behavior, with equal weight for
  each target-actor/session stratum.
- Target-behavior comparisons separately for self and other. Target labels
  are permuted only within the same target actor and session.
- Post-offset versus pre-onset baseline tests separately for self and other.
  Sign flips are clustered by target actor and session.
- Window-level and bin-wise CSV tables. Bin-wise tables include max-T
  family-wise corrected p-values.

The analysis output and figures are written only to the separate actor-pipeline
folder above.

## Check

`test_actor_behavior_probability.m` uses a synthetic four-animal session. It
checks target assignment, empty actor blocks, self behavior, the logical OR
across multiple other animals, summaries, paired-target eligibility, and
baseline result dimensions.
