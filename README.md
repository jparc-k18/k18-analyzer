# k18analyzer-e05-shared

Analyzer framework for the J-PARC E05 experiment (Oct-Nov 2015,
K1.8/SKS), based on the K18 analyzer v20150402.

This is the **shared** distribution: the analyzer code and the
v20150402 baseline parameter set. The calibration constants and
the run bookkeeping derived from the E05 data are not included
-- see "What this repository does not contain" below.

## Build

```
make            # lib + usr + dst -> bin/
```

The HD DAQ unpacker is an external dependency. It defaults to
`$(HOME)/work/k18analyzer/e05/unpacker`; override it if yours
lives elsewhere:

```
make unpacker_install=/path/to/unpacker
```

## Setup

The following are machine-specific and are not tracked. Create
them as symlinks (or plain directories) in the repository root:

| name         | points at                                  |
|--------------|--------------------------------------------|
| `e05_oct`    | the E05 Oct-2015 raw data directory        |
| `root`       | where the output ROOT files should go      |
| `fieldmap`   | the SKS field map directory                |
| `K18LogBook` | the K1.8 logbook directory (optional)      |

## Layout

- `src/`, `include/` -- analyzer library (tracking, hodoscope,
  Runge-Kutta, parameter managers).
- `usr/` -- event-loop analyzers (`bin/<Name>` reads a conf file
  and a raw data file). Every analyzer writes `runnum`, `evnum`
  and `trigflag` and fills its tree for EVERY event, so all trees
  of one run have identical event counts and can be matched by
  `evnum`.
- `dst/` -- post-stage programs combining analyzer trees.
  `DstSksHodoscope` matches SksTracking x SksHodoscope by `evnum`
  and RECOMPUTES the TOF times/dE from the raw TDC/ADC with the
  current HodoParam/HodoPHC, so the TOF calibration iterates
  without re-running the Runge-Kutta tracking. `DstKKAna` /
  `DstKPiAna` build the (K-,K+) / (K-,pi+) pair trees and the
  missing-mass histograms.
- `param/` -- the v20150402 baseline parameter set, plus the E05
  channel maps, digit maps and unpacker configuration.
- `runmanager/` -- LSF batch manager. It is driven by run lists
  (yaml); `runmanager/runlist/example.yml` documents the format.

## What this repository does not contain

The analyzer runs on a per-run conf file that points at per-run
calibration tables:

```
param/conf/e05/analyzer_<run>.conf
param/HDPRM/e05/HodoParam_<run>      TDC/ADC conversion
param/HDPHC/e05/HodoPHC_<run>        slewing correction
param/DCTDC/e05/DCTdcCalib_<run>     drift-chamber TDC offsets
param/DCDRFT/e05/DCDriftParam_<run>  drift time-to-length
param/DCGEO/e05/DCGeomParam_e05      DC alignment
param/FIELD/e05/field_table.txt      per-run field records
param/MOMCALIB/e05/mom_calib.json    momentum calibration
```

**None of these are distributed here.** They are the output of
the E05 calibration chain and are held separately. `param/*/e05/`
is in `.gitignore`, so dropping a calibration set into place does
not dirty the working tree.

To run without them, write your own conf modelled on
`param/conf/analyzer_20151025.conf`, which references only
baseline tables that ARE included.

The E05 run lists are not distributed either: which runs belong
to which target, beam momentum and trigger setting is bookkeeping
that goes with the calibration. `runmanager/runlist/` therefore
holds only `example.yml`, and `*.yml` and subdirectories under it
are gitignored, so your own run lists can live there without
dirtying the working tree.

Related consequences:

- `dst/DstKKAna.cc` applies a momentum calibration through the
  `H_COEF` / `K_COEF` / `PO0` / `PO1` constants. They are **zeroed
  here**, so the momenta pass through unchanged. Fit them for your
  own data set before trusting the missing-mass spectra.
- The calibration macros and the step-by-step calibration guides
  are not part of this repository.

## Run lists

```sh
python3 runmanager/run.py runmanager/runlist/my-runlist.yml
```

Every path in a run list is relative to the repository root, so a
run list needs no absolute path anywhere:

```yaml
DEFAULT:
  queue: l      # use "l", not "s", for mass production
  bin:   bin/Hodoscope
  conf:  param/conf/my.conf
  data:  e05_oct
  root:  root/v1
RUN:
  10217:        # reads e05_oct/run10217.dat.gz
```

An optional `WORKDIR` key overrides that root; set it only if the
run list drives an analyzer outside this checkout. See
`runmanager/runlist/example.yml` for the full key reference.
