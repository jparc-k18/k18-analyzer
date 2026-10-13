# Changes from k18analyzer v20150402

Baseline: `k18analyzer_20150402.tar.gz`
(`/group/had/sks/E07/software/`), the K1.8 offline analyzer
snapshot this repository started from. That tarball ships code
only -- `Makefile`, `include/`, `src/`, `usr/`, `sh/`, `conf/` --
so everything under `param/`, `runmanager/` and `dst/` here is an
addition rather than a modification.

## Summary

| area | v20150402 | here | note |
|------|-----------|------|------|
| `include/*.hh` | 76 | 66 | 10 removed, 10 modified, 56 as-is |
| `src/*.cc` | 59 | 49 | 10 removed, 14 modified, 35 as-is |
| `usr/*.cc` | 13 | 11 | 2 removed, 10 modified, 1 as-is |
| `src/backup`, `usr/backup` | 15 | -- | dropped |
| `sh/` | 14 | -- | dropped (E13-era job scripts) |
| `conf/` | 7 | -- | dropped (E13 conf files, not relocated) |
| `dst/` | -- | 4 | new post-stage programs |
| `param/` | -- | 635 | baseline parameter set, added |
| `runmanager/` | -- | 17 | LSF batch manager, added |

## 1. CLHEP dropped in favour of ROOT

`ThreeVector` / `LorentzVector` were `CLHEP::Hep3Vector` and
`CLHEP::HepLorentzVector`; they are now `TVector3` and
`TLorentzVector`, and `clhep-config` is gone from the `Makefile`.
The analyzer therefore builds against ROOT alone.

The migration is mostly mechanical (`.mag()` -> `.Mag()`,
`setX()` -> `SetX()`) across `Kinematics.cc`, `DCAnalyzer.cc`,
`RungeKuttaUtilities.cc`, `SksTrack.cc`, `SksParticle.cc`,
`FieldMan.cc`, `DopplerCorr.cc` and the `usr/` analyzers, with two
semantic traps handled explicitly:

- `TVector3` has no `operator/(double)` and no `operator<<`;
  both are added as free functions in `ThreeVector.hh`.
- `Hep3Vector::rotateZ()` returns a rotated copy, while
  `TVector3::RotateZ()` rotates in place and returns `void`.
  The four call sites in `EvDisp.cc` are rewritten accordingly --
  the naive substitution silently yields the unrotated vector.

## 2. Hyperball-J and SksMinus detectors removed

E05 has no germanium or PWO array, so the Hyperball-J readout was
taken out rather than left dead: `GeRawHit`, `GeHit`,
`GeAdcCalibMan`, `GeTdcMan`, `PWORawHit`, `PWOHit`, `PWOTdcMan`,
`HBJHit`, `HBJAnalyzer` and the SksMinus `SP0Hit` (10 headers +
10 sources), together with `usr/UserHyperballJ.cc` and
`usr/UserKPiGAna.cc`.

Consequently `DetectorID.hh`, `RawData.{hh,cc}`,
`HodoAnalyzer.{hh,cc}`, `SksLib.hh` and `ConfMan.{hh,cc}` lose
their Ge/PWO/SP0 members, and the `GeADC:`, `GeTDC:`, `PWOTDC:`
conf keys are no longer recognised.

## 3. Bug fixes in the inherited code

These are defects in v20150402 itself, not E05 adaptations.

**Runge-Kutta (`RungeKuttaUtilities`)** -- three transcription
errors in the transport-derivative code:

| location | was | is |
|----------|-----|----|
| `RKcalcDeltaFieldIntegral` | `dkzq = intg.ky/...` | `intg.kz` |
| `RKtraceToLast` | `hpini.dUdZ()` passed as `dVdZ` | `hpini.dVdZ()` |
| `RKDeltaFieldIntegral` ctor | `dudq(Dydq)` | `dudq(Dudq)` |

Each puts a wrong element into the derivative matrix, so the
covariance/fit convergence is degraded rather than obviously
broken -- which is why they survive silently.

**Missing return values** -- `DCAnalyzer::DecodeRawHits()`,
`DCAnalyzer::TDCrearrangeSDCout()`, `ConfMan::InitializeEvDisp()`
and the `const char*` / `std::string` overloads of
`DCGeomMan::Initialize()` and `FieldMan::Initialize()` fell off
the end of a non-`void` function (the overloads also discarded
the result of the real `Initialize()`).

**Uninitialised members** -- `ConfMan::ConfMan(const std::string&)`
initialised only the two file-name members, leaving every manager
pointer and the momentum/field scalars indeterminate; the
default constructor did initialise them. It now has the same
initialiser list.

**BcOut -> SKS frame link** -- `DCLocalTrack::DoFitBcSdc()` and
`LocalTrackSearchBcOutSdcIn()` placed the K18 target at
`GetLocalZ(130)` alone, ignoring `GetLocalZ(0)` (-145.3 mm), the
target position in the SKS frame. The whole BcOut arm was
displaced by that amount relative to SdcIn. Both now use
`GetLocalZ(130) - GetLocalZ(0)`.

## 4. Tracking changes for E05

**BcOut x SdcIn combination search** -- the original
`LocalTrackSearchBcOutSdcIn()` required *exactly one* cluster on
*every* pair plane and returned zero otherwise. With BC3 layers
1-2 dead in E05 that rejects every event. It now enumerates
combinations with `makeindex()` like the standard
`LocalTrackSearch()`, guarded by a 1e5 cap on the combination
count, and the minimum layer requirement is 16 of 20 live layers
(was 19).

**SKS initial momentum** -- the `ver1_Sks` path seeded the
Runge-Kutta fit with a fixed momentum, which is up to 50% off
across the E05 settings (1.1-1.9 GeV/c). The seed is now the
empirical dispersion `1/p = c0 + c1*uIn + c2*uOut` fitted over
the beam-through scan and the CH2 production tracks, good to
about +-3% and independent of the run condition.

## 5. Output tree conventions

Every `usr/` analyzer now writes `runnum` and `evnum` branches
(from `GUnpacker`) and **fills its tree for every event**: each
early `return true` in the event loop became
`{ tree->Fill(); return true; }`. All trees of one run therefore
have identical entry counts and can be matched by `evnum` -- this
is what makes the `dst/` stage possible.

Two things had to change to make that safe:

- Size branches are initialised to `0`, not `-1`.
  `TTree::Fill()` crashes on a negative array size when the fill
  section is skipped.
- In the six tracking analyzers, track and hit multiplicities
  are clamped to `MaxHits` before being stored, so a noisy event
  cannot overrun the fixed-size arrays.

**Trigger flag off-by-one** -- `event.trigflag[seg]` was written
with the 1-based segment number, so `trigflag[0]` was never
filled and the last trigger bit was dropped. It is now
`trigflag[seg-1]`. `DetectorID.hh` gains a `TrigFlag` enum naming
all 20 bits, decoded from K18LogBook41 p.199 and cross-checked
against the per-bit firing rates in the data.

**Scattered-particle timing** (`UserSksTracking`) -- `m2` used the
hardcoded `offset = 1.3108917166` from `Parameter.hh`, an E13-era
constant for kaons at 1.8 GeV/c. The Runge-Kutta path length
starts at the target, so the BH2 -> target beam leg is now
subtracted explicitly, from the surveyed `DCGEO` positions
(`K18Target` id 130, `BH2` id 202) and the beam velocity implied
by `PK18:`. `UserKPiAna` and `UserSksHodoscope` still use the old
`offset` (the latter reads it from the `USER:` parameter file);
`Parameter.hh` itself is untouched.

**Hodoscope calibration histograms** (`UserHodoscope`) -- added
per-segment BTOF/CBTOF vs dE 2-D histograms for BH1 and BH2, the
input to the slewing-correction fits.

## 6. Reproducibility metadata

`ConfMan` keeps every two-token `KEY VALUE` line of the conf file
(not just the whitelisted keys) and `ConfMan::WriteParameters()`,
called from `src/Main.cc`, writes into the output ROOT file:

- a top-level `git` `TNamed` holding `git log -1`,
- a `conffile` `TNamed`,
- a `param` directory with one `TMacro` per referenced parameter
  file (the file contents, embedded) and one `TNamed` per scalar.

So a ROOT file documents the exact conditions that produced it.
Two guards matter:

- files larger than 1 MB are recorded by path only. The SKS field
  map is ~72 MB and `hadd` cannot merge or dedupe `TMacro`
  objects, so embedding it made a merged output balloon to
  several GB of duplicated copies.
- only the `skip == 0` sub-job writes the metadata, otherwise a
  unit-split run accumulates one copy per chunk.

## 7. New `dst/` stage

Post-stage programs that combine Stage-1 trees by `evnum`:

- `DstSksHodoscope` -- matches `SksTracking` x `SksHodoscope` and
  **recomputes** the TOF times and dE from the raw TDC/ADC with
  the current `HodoParam`/`HodoPHC`, so the TOF calibration can
  iterate without re-running the Runge-Kutta tracking.
- `DstKKAna` / `DstKPiAna` -- (K-,K+) / (K-,pi+) pair trees:
  vertex, closest distance, lab angle, `m2` and the missing-mass
  histograms. No physics cut is applied to the tree; the
  selection stays visible at the plotting stage. The target
  material for the energy loss comes from a `TARGET:` conf key.
- `DstMetaUtil.hh` -- the `ConfMan::WriteParameters()` equivalent
  for the Dst stage, which additionally copies the metadata
  already embedded in each Stage-1 input.

## 8. Build system

- CLHEP removed (see 1).
- `dst` added as a third build target next to `lib` and `usr`.
- `unpacker_install` is `?=`, so the unpacker location can be
  overridden from the environment or the command line.
- `lib/` and `bin/` are created on demand; both are gitignored,
  so without this a fresh clone fails at the archive step.
- The build file is tracked as `Makefile.org` and `Makefile` is
  gitignored, matching the other branches: `cp Makefile.org
  Makefile` before the first build.

## 9. Conf file keys

Two spellings were in circulation: the v20150402 code reads
`DCGEO:` / `DCTDC:` / `DCDRFT:`, while the E05-era parameter
directory writes `DCGEOM:` / `TDCCALIB:` / `DRFTPM:`. The
canonical keys here follow the **current K1.8 analyzer**
(HYPS/E70), which keeps the short forms and shortens
`UNPACKER:` to `UNPACK:`:

| canonical (here, HYPS/E70) | v20150402 | E05 param dir |
|----------------------------|-----------|---------------|
| `UNPACK:`  | `UNPACKER:` | `UNPACKER:` |
| `DCGEO:`   | `DCGEO:`    | `DCGEOM:`   |
| `DCTDC:`   | `DCTDC:`    | `TDCCALIB:` |
| `DCDRFT:`  | `DCDRFT:`   | `DRFTPM:`   |

Every other key (`DIGIT:`, `CMAP:`, `HDPRM:`, `HDPHC:`, `USER:`,
`K18TM:`, `PK18:`, `FLDMAP:`, `FLDNMR:`, `FLDCALC:`) already
matches.

The three right-hand spellings are accepted as aliases and
normalised on read, so an old conf file still works and the key
recorded in the output ROOT file is the canonical one either way.
`runmanager` also reads the unpacker path out of the conf file
and now looks for `UNPACK`, the name it uses upstream.

`GeADC:`, `GeTDC:` and `PWOTDC:` are no longer recognised (see 2)
and `TARGET:` is new (see 7).

One inherited oddity is left alone: `USER:` appears in the E05
conf files but nothing in this lineage consumes it -- the
`FILEMAN:` key is what feeds `FileManager`. In the current K1.8
analyzer `USER:` drives `UserParamMan`, which does not exist
here.

## 10. Not a code change

`param/` (the K1.8 parameter directory: the inherited baseline
tables plus the E05 channel and digit maps) and `runmanager/`
(Y. Nakada's LSF batch manager, v4.1) come from outside the
tarball and are carried here essentially unmodified -- the
exceptions are the work-directory default described in the git
history, and three small fixes to inherited parameter files
(a malformed `737.5.0` field in `HodoParam_e05_0`, the
`disagreement local` flag in `unpacker.xml`, and the `../`
relative paths in `analyzer_20151025.conf`).

The E05 calibration output and run bookkeeping are not part of
this branch -- see README.md.
