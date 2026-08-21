# Runlist for E05 Oct 2015

Runlist yaml files for `runmanager/run.py`.
See `example.yml` for the yaml format.

## hodo/

Hodoscope runlists per run set, generated from
`e05_oct/misc/comment.txt`.
Runs without a data file in `e05_oct/` are excluded
(10011-10017, 10058-10061, 10134, 10239-10269, 10328, 10484-10485).
Runs with fewer than 1000 events in `recorder.log` are also excluded
as empty. Dummy runs, empty-comment runs and beam-tuning junk runs
are not assigned to any set.

Usage:

```sh
python3 runmanager/run.py runmanager/runlist/hodo/hodo-prod-c-m18-2.yml
```

### Production run sets

| File | Description | Runs |
|------|-------------|------|
| hodo-prod-ch2-p18.yml | CH2 target @+1.8 GeV/c | 10050-10056 |
| hodo-prod-ch2-p17.yml | CH2 target @+1.7 GeV/c | 10069-10083 |
| hodo-prod-ch2-p16.yml | CH2 target @+1.6 GeV/c | 10093-10104 |
| hodo-prod-ch2-p19.yml | CH2 target @+1.9 GeV/c | 10107-10118 |
| hodo-prod-ch2-p15.yml | CH2 target @+1.5 GeV/c | 10120-10130 |
| hodo-prod-ch2-m18.yml | CH2 target @-1.8 GeV/c | 10138-10176 |
| hodo-prod-ch2thick-m18.yml | CH2 4.8 g/cm2 @-1.8 GeV/c | 10178-10182 |
| hodo-prod-c-m18-1.yml | C target @-1.8 GeV/c (KinxTOF 1/150) | 10217-10322 |
| hodo-prod-c-m18-2.yml | C target @-1.8 GeV/c (KinxTOF 1/140) | 10346-10605 |
| hodo-sigma-p15.yml | Sigma- run @+1.5 GeV/c (SKS 400 A) | 10143-10145 |
| hodo-sigma1385-m15.yml | Sigma(1385) w/ C @-1.5 GeV/c | 10545-10548 |
| hodo-pik-study-m18.yml | (pi-, K+) study w/ C @-1.8 GeV/c | 10550-10552 |
| hodo-empty-target-m18.yml | Empty target @-1.8 GeV/c | 10522-10527 |

### Beam through / profile check sets

| File | Description | Runs |
|------|-------------|------|
| hodo-beamthrough-1101.yml | Beam through (CH2 / empty / 12C, 11/01) | 10186-10215 |
| hodo-beamthrough-1115.yml | Beam through w/o SKS field (11/15) | 10492-10519 |
| hodo-profile-1115.yml | Beam profile check (11/15) | 10490-10521 |
| hodo-counter-check.yml | Profile and counter check (10/27) | 10064-10067 |
| hodo-bc3-check.yml | BC3 noise and profile check | 10132, 10483, 10486 |

### Tuning / calibration sets

| File | Description | Runs |
|------|-------------|------|
| hodo-pedestal.yml | Pedestal runs | 10009-10606 (12 runs) |
| hodo-ktune-p18.yml | K tuning @+1.8 GeV/c | 10018-10033 |
| hodo-ktune-p16.yml | K tuning @+1.6 GeV/c | 10089-10092 |
| hodo-ktune-m18.yml | K tuning @-1.8 GeV/c | 10137 |
| hodo-ac-study.yml | AC threshold study | 10034-10042 |
| hodo-tof-study.yml | TOF Vth / efficiency study | 10043-10048 |
| hodo-scaler-study.yml | Scaler study / check | 10495-10533 |
| hodo-kic-test.yml | KIC test and cosmic ray | 10324-10334 |
| hodo-pwo-cosmic.yml | PWO study and cosmic ray | 10327-10487 |
