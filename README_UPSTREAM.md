# K1.8 native Runge-Kutta momentum reconstruction

## Summary

This tree adds the modules needed to reconstruct K1.8 momentum with a
native-coordinate Runge-Kutta transport through QQDQQ, to reconstruct K18 and
S-2S tracks directly from a `g4s2s` tree, and to evaluate generator-6381
missing mass from the entry-aligned tracking outputs. The native K18 fit
propagates a VO track backward and varies its momentum until the calculated
BFT coordinate matches the measured cluster position.

The established D2U transport-matrix result remains the primary K1.8 momentum
for standalone resolution and missing-mass results.  The native RK result is a
field-model diagnostic and cross-check.

## Configuration

- `K18FLDMAP` selects the regular QQDQQ field map.
- `K18FLDNMR/K18FLDCALC` applies the common map scale.
- `K18RKCharge` sets the transported particle charge sign.
- `K18NativePMin/PMax` defines the fit interval.
- `K18NativeRKStep` and `K18NativeRKDriftStep` set field and drift step sizes.
- `K18NativeFitMaxIteration`, `K18NativeFitMomentumTolerance`, and
  `K18NativeFitResidualTolerance` control convergence.
- `K18NativeMaxPath` and `K18NativeBFTResidualMax` reject invalid transport or
  poor BFT matches.
- `K18GlobalScale` and individual magnet scale keys apply when the analytic
  field model is used instead of a map.
- `S2sGeant4AcceptedPdg` and `S2sGeant4AcceptSecondaries` select the G4 track
  used for S-2S hit reconstruction.
- `S2sRKPolarity` sets the S-2S RK charge sign.
- `DstS2sTrackingSeedMode=1` uses `PK18`; mode 2 uses the fixed
  `DstS2sTrackingSeedMomentum`. Neither mode reads eventwise truth momentum.

Input track positions are millimetres and momentum is GeV/c in the K1.8 native
frame.  Field-map file coordinates are centimetres and field values are Tesla.

## File guide

| File | Purpose |
| --- | --- |
| `include/K18BeamlineRK.hh` | Defines propagated states, fit results/status, configuration accessors, and the RK fitting interface. |
| `src/K18BeamlineRK.cc` | Implements QQDQQ coordinates and fields, RK4 propagation, BFT crossing interpolation, bracketed momentum solving, bounded minimization, and track-result filling. |
| `include/K18FieldMap.hh` | Defines a validated regular-grid field-map interface and map metadata accessors. |
| `src/K18FieldMap.cc` | Loads every grid point with duplicate/missing-record checks and evaluates the field by trilinear interpolation. |
| `dst/DstK18TrackingGeant4.cc` | Reconstructs K18 D2U/TM momentum and the optional native-RK diagnostic directly from `g4s2s`. |
| `dst/DstS2sTrackingGeant4.cc` | Reconstructs S-2S local and RK tracks from G4 detector hits with deterministic smearing and explicit particle/polarity selection. |
| `dst/DstMissingMassGeant4.cc` | Reads entry-aligned `g4s2s`, `k18track`, and `s2s`; writes the `bbflow` tree, missing-mass/component histograms, and cutflow/fit CSV. |
| `runmanager/runlist/README_K18_GEANT4.md` | Gives 6380/6381 analyzer configuration templates and commands for the three DST stages. |

## Geant4 DST chain

Build and usage commands are in
[`runmanager/runlist/README_K18_GEANT4.md`](runmanager/runlist/README_K18_GEANT4.md).
Generator 6380 is used for beam-through tracking validation. Generator 6381
adds reaction truth and can be processed by `DstMissingMassGeant4`. The
missing-mass DST uses D2U/TM `p_3rd` as the primary K18 momentum and supports
configurable beam/scattered masses and PDGs, slab/cylinder target geometry,
target/cell/shell Bethe--Bloch parameters, selection cuts, and output paths.

## Verification

- `make -f Makefile.org -j2 dst` succeeds.
- On the 100000-entry generator-6381 K- to pi- sample, current pro and E63
  tracking outputs have identical K18/S-2S entries, track multiplicities, and
  reconstructed momenta.
- `DstMissingMassGeant4` selected 268 smoke-test events with zero event-number
  mismatches. Its CSV is byte-identical to `AnaReactionBBFlowBudget`, and the
  checked missing-mass and vertex histograms have zero bin-content difference.
  This verifies reproduction only; the selected count is not a final
  high-statistics physics result.
