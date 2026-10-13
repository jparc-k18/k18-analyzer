---
stylesheet: https://cdnjs.cloudflare.com/ajax/libs/github-markdown-css/2.10.0/github-markdown.min.css
body_class: markdown-body
pdf_options:
  format: A4
  margin: 24mm 16mm
  displayHeaderFooter: true
  headerTemplate: |-
    <style>
      section {
        margin: 0 auto;
        font-family: system-ui;
        font-size: 11px;
      }
    </style>
    <section>
    </section>
  footerTemplate: |-
    <section>
      <div>
        <span class="pageNumber"></span>
        <!-- /<span class="totalPages"></span> -->
      </div>
    </section>
---
K1.8 analyzer README (E05)
==========================

<div style="text-align: right;">
Rev. 2026.10.13<br>
</div><br>

Analyzer for the J-PARC E05 experiment (Oct-Nov 2015, K1.8/SKS),
based on the K18 analyzer v20150402.

It is assumed to work on KEKCC, Red Hat Enterprise Linux release
9.3 (Plow). Use gcc 11.4.1 and root 6.32.04.

This branch does NOT reproduce the results of v20150402: the
tracking and the output tree conventions have changed. See
diff.md for the full list.

## Environment setting

See KEKCC.pdf and set environment variables.
Anaconda/Python environment is needed to use runmanager.

## Unpacker

The unpacker compiled with gcc 11.4.1 is placed in the group
directory, /group/had/sks/software/unpacker/unpacker.jparc-k18.
This is updated constantly.

If you want to install in local, install as follows.

```sh
$ git clone https://github.com/jparc-k18/unpacker.git
$ cd unpacker/src
$ cp Makefile.org Makefile
$ make
```

Check if the "unpacker-config" command is available.

```sh
$ unpacker-config --version
```

## K1.8 analyzer

Install the K1.8 analyzer.

```sh
$ git clone https://github.com/jparc-k18/k18-analyzer.git
$ cd k18-analyzer
$ git checkout e05
$ cp Makefile.org Makefile
$ make
```

The Makefile picks up the group unpacker by default. Override it
if yours lives elsewhere.

```sh
$ make unpacker_install=/path/to/unpacker
```

e.g.) Hodoscope,
Usage: Hodoscope [analyzer config file] [data input stream]
[output root file]

```sh
$ ./bin/Hodoscope param/conf/analyzer_20151025.conf \
/hsm/had/sks/E05/SKS2015Oct/e05_oct/run10217.dat.gz hoge.root
```

## Setup

The following are machine-specific and are not tracked. Create
them as symlinks (or plain directories) in the repository root.

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
- `param/` -- the E05-era baseline parameter set. Unlike the
  recent branches, this one carries `param/` in tree: the E05 set
  predates the current jparc-k18/param layout and includes
  directories that have since been pruned there.
- `runmanager/` -- LSF batch manager, driven by run lists (yaml);
  `runmanager/runlist/example.yml` documents the format.

## What this branch does not contain

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

__None of these are distributed here.__ They are the output of
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
  `H_COEF` / `K_COEF` / `PO0` / `PO1` constants. They are
  __zeroed here__, so the momenta pass through unchanged. Fit
  them for your own data set before trusting the missing-mass
  spectra.
- The calibration macros and the step-by-step calibration guides
  are not part of this branch.

### runmanager

__runmanager__ is a script for managing jobs on KEKCC.
Prepare runlist.yml by referring to the example.yml.
Note that the indentation determines the nest depth in yaml.

Every path in a run list is relative to the repository root, so a
run list needs no absolute path anywhere:

```yml
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

Run run.py.
Note that the process runs until all jobs have finished.
Ctrl-C kills all jobs and terminates the process.

```sh
$ ./runmanager/run.py runmanager/runlist/foo.yml
```

Run monitor.py on another tty to see the progress of the jobs.
The job status is updated in the "stat" directory, using the same
name as the runlist in json format.

```sh
$ ./runmanager/monitor.py runmanager/stat/foo.json
```

![runmanager.png](runmanager.png "caption")
