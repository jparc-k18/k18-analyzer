# K18 Geant4 tracking and missing-mass examples

## Reaction configuration set

The five generator-6381 reactions below have matching Geant4, K18, and S-2S
configurations. Run Geant4 from the `k18-geant4` repository and both tracking
stages from this repository.

| Tag | Reaction | Geant4 configuration | Analyzer configurations |
| --- | --- | --- | --- |
| `e90_cusp` | `d(K-,pi-)Sigma-N` cusp | `k18_6381_e90_cusp_example.conf` | `k18_6381_e90_cusp_{k18,s2s}.conf` |
| `e75_phase1_li7` | `7Li(K-,K+)7_XiH` | `k18_6381_e75_phase1_li7_example.conf` | `k18_6381_e75_phase1_li7_{k18,s2s}.conf` |
| `kpi_c12lambda` | `12C(K-,pi-)12_LambdaC` | `k18_6381_kpi_c12lambda_example.conf` | `k18_6381_kpi_c12lambda_{k18,s2s}.conf` |
| `kk_c12xibe` | `12C(K-,K+)12_XiBe` | `k18_6381_kk_c12xibe_example.conf` | `k18_6381_kk_c12xibe_{k18,s2s}.conf` |
| `pik_c12lambda` | `12C(pi+,K+)12_LambdaC` | `k18_6381_pik_c12lambda_example.conf` | `k18_6381_pik_c12lambda_{k18,s2s}.conf` |

Before use, replace the `/path/to/...` entries in the analyzer configurations
with an equivalent analyzer parameter release and the local `k18-geant4`
checkout. S-2S field maps are read directly from
`/group/had/sks/fieldmap/S2S` on KEKCC. The older `k18_638*_example.conf`
files use the same placeholder convention.

## Build and run

Build the three stages with:

```sh
make -f Makefile.org bin/DstK18TrackingGeant4 \
  bin/DstS2sTrackingGeant4 bin/DstMissingMassGeant4
```

For example, run K18 and S-2S tracking on the same E90 generator ROOT file:

```sh
bin/DstK18TrackingGeant4 \
  runmanager/runlist/k18_6381_e90_cusp_k18.conf \
  e90_g4.root e90_k18.root
bin/DstS2sTrackingGeant4 \
  runmanager/runlist/k18_6381_e90_cusp_s2s.conf \
  e90_g4.root e90_s2s.root
```

Always keep G4, K18, S-2S, and analysis outputs on distinct canonical paths.
Before analysis, run the E63 stage validator:

```sh
python /path/to/k18-analyzer-e63/scripts/validate_mc_stage_contract.py \
  --g4 e90_g4.root --s2s e90_s2s.root
```

Generator 6380 is a beam-through sample and stops at tracking validation.
Generator 6381 also provides `ReactionBeamVertexTransport`,
`ReactionBeamVertex`, and `ReactionScat`, allowing the reconstructed missing
mass to be evaluated directly:

## `DstMissingMassGeant4` parameter presets

Use the following reaction-specific values with the common command below.
Masses are in GeV/c2, lengths in mm, density in g/cm3, and mean excitation
energy `I` in eV.

| Tag | beam/scat PDG | beam/scat mass | target mass | target geometry and material |
| --- | --- | --- | ---: | --- |
| `e90_cusp` | `-321/-211` | `0.493677/0.13957039` | `1.87561294257` | cylinder `R=27`, half-`Y=50`; `rho=0.169`, `I=19.2`, `Z/A=0.496524329692`; plane `Z=-5114` |
| `e75_phase1_li7` | `-321/321` | `0.493677/0.493677` | `6.533832826901` | slab `16.853932584`, `100x100`; `rho=0.534`, `I=40.0`, `Z/A=0.432214378331`; plane `Z=-5023.5` |
| `kpi_c12lambda` | `-321/-211` | `0.493677/0.13957039` | `11.174863235340` | C slab `5.555555556`, `200x100`; `rho=1.8`, `I=81.0`, `Z/A=0.499542086420`; plane `Z=-5023.5` |
| `kk_c12xibe` | `-321/321` | `0.493677/0.493677` | `11.174863235340` | C slab `52.0`, `200x100`; otherwise as above |
| `pik_c12lambda` | `211/321` | `0.13957039/0.493677` | `11.174863235340` | C slab `5.555555556`, `200x100`; otherwise as above |

```sh
bin/DstMissingMassGeant4 \
  --g4=TAG_g4.root --k18=TAG_k18.root --s2s=TAG_s2s.root \
  --root=TAG_missing.root --csv=TAG_missing.csv --label=TAG \
  --beam-pdg=BEAM_PDG --scat-pdg=SCAT_PDG \
  --beam-mass=BEAM_MASS --scat-mass=SCAT_MASS \
  --target-mass=TARGET_MASS --target-shape=slab \
  --target-thickness-mm=THICKNESS --target-density=DENSITY \
  --target-i-ev=I --target-z-over-a=Z_OVER_A \
  --target-size-x-mm=SIZE_X --target-size-y-mm=SIZE_Y \
  --target-plane-z-mm=TARGET_Z --target-z-margin-mm=0 \
  --theta-min-deg=2 --theta-max-deg=14
```

For E90, replace the slab geometry options with:

```text
--target-shape=cylinder --target-radius-mm=27
--target-half-length-y-mm=50 --target-thickness-mm=54
--target-size-x-mm=54 --target-size-y-mm=100
--target-plane-z-mm=-5114 --target-z-margin-mm=0
```

The validated E90 material budget additionally used a 0.25-mm Mylar target
cell and 0.25-mm G10 vacuum shell; pass the corresponding `--cell-*` and
`--vacuum-shell-*` options when reproducing that budget.

`DstMissingMassGeant4` keeps the `AnaReactionBBFlowBudget` output contract:
the branch-limited `bbflow` tree, component histograms, fit/cutflow CSV, and
the reconstructed-vertex plus mean Bethe--Bloch correction flow. Its primary
K18 momentum is D2U/TM `p_3rd`; native K18 RK remains diagnostic.

The DST rejects canonical path collisions, missing required branches, entry
count mismatches, and returns status 2 for event-number mismatches. Also inspect
underflow and overflow before quoting a core resolution.
