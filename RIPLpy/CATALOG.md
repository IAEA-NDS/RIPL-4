# RIPLpy Data Catalog

Machine-generated schema for every RIPLpy database (regenerate with
`python tools/generate_schema.py`). The authoritative machine-readable
form is `SCHEMA.json`; this file is the human-readable rendering.

## densities

### bsfg
_legacy/absent in this release_

### bsk14_comb  (8508 entries, key=Nuclide, entry=Entry)

| Field | Type | Description |
|-------|------|-------------|
| `n` | Nuclide | Target nucleus |
| `positive_parity` | riplpy.densities.hfb.ParityData | Level density data for positive parity states |
| `negative_parity` | riplpy.densities.hfb.ParityData | Level density data for negative parity states |

### bskg3_comb  (7677 entries, key=Nuclide, entry=Entry)

| Field | Type | Description |
|-------|------|-------------|
| `n` | Nuclide | Target nucleus |
| `positive_parity` | riplpy.densities.hfb.ParityData | Level density data for positive parity states |
| `negative_parity` | riplpy.densities.hfb.ParityData | Level density data for negative parity states |

### ct
_legacy/absent in this release_

### egsm  (291 entries, key=Nuclide, entry=Entry)

| Field | Type | Description |
|-------|------|-------------|
| `n` | Nuclide | Target nucleus |
| `Io` | float | Ground state spin of target nucleus |
| `Bn` | float | Neutron binding energy / Qn [MeV] |
| `Do` | float | Experimental s-wave resonance spacing Dobs [keV] |
| `Derr` | float | Uncertainty on resonance spacing dDobs [keV] |
| `Esh` | float | Shell correction energy [MeV] |
| `Def` | float | Ground-state deformation parameter |
| `Dcalc` | float | Calculated s-wave resonance spacing [keV] |
| `dap` | float | Upper uncertainty on level density parameter [MeV^-1] |
| `a` | float | Level density parameter at Bn [MeV^-1] |
| `dam` | float | Lower uncertainty on level density parameter [MeV^-1] |
| `a_sys` | float | Systematic level density parameter a_sys [MeV^-1] |
| `a_ratio` | float | Ratio a_exp / a_sys |

### egsm_norm  (82 entries, key=Element, entry=Entry)

| Field | Type | Description |
|-------|------|-------------|
| `element` | Element | Target element |
| `factor` | float | Normalization factor (a_exp/a_sys averaged over isotopes) |

### hfb  (0 entries, key=Nuclide, entry=None)

_(no field metadata)_

### qrpabe  (2706 entries, key=Nuclide, entry=Entry)

| Field | Type | Description |
|-------|------|-------------|
| `n` | Nuclide | Target nucleus |
| `positive_parity` | riplpy.densities.hfb.ParityData | Level density data for positive parity states |
| `negative_parity` | riplpy.densities.hfb.ParityData | Level density data for negative parity states |

### shellcorr_mk
_legacy/absent in this release_

### shellcorr_ms  (6830 entries, key=Nuclide, entry=ExtendedEntry)

| Field | Type | Description |
|-------|------|-------------|
| `n` | Nuclide | Target nucleus |
| `shell` | float | Shell correction energy [MeV] |
| `corr` | float | Deformation correction energy [MeV] |
| `beta2` | float | Quadrupole deformation parameter |
| `beta4` | float | Hexadecapole deformation parameter |

### thfb_comb  (7334 entries, key=Nuclide, entry=Entry)

| Field | Type | Description |
|-------|------|-------------|
| `n` | Nuclide | Target nucleus |
| `positive_parity` | riplpy.densities.hfb.ParityData | Level density data for positive parity states |
| `negative_parity` | riplpy.densities.hfb.ParityData | Level density data for negative parity states |

## fission

### bskg3_barriers  (2449 entries, key=Nuclide, entry=Entry)

| Field | Type | Description |
|-------|------|-------------|
| `n` | Nuclide | Target nucleus |
| `inner` | dict | Inner barrier {E[MeV], B20, B22, B30} |
| `outer1` | dict | First outer barrier {E[MeV], B20, B22, B30} |
| `outer2` | dict | Second outer barrier {E[MeV], B20, B22, B30} |
| `isomer` | dict | First shape isomer {E[MeV], B20, B22, B30} |
| `isomer2` | dict | Second shape isomer {E[MeV], B20, B22, B30} (None on the legacy layout) |

### d1m_barriers  (45 entries, key=Nuclide, entry=Entry)

| Field | Type | Description |
|-------|------|-------------|
| `n` | Nuclide | Target nucleus |
| `inner` | dict | Inner barrier {E[MeV], B20, B22, B30} |
| `outer` | dict | Outer barrier {E[MeV], B20, B22, B30} |

### empire_barriers  (76 entries, key=Nuclide, entry=Entry)

| Field | Type | Description |
|-------|------|-------------|
| `n` | Nuclide | Target nucleus |
| `symmetry` | str | Saddle point symmetry designation |
| `Va` | float | First barrier height [MeV] |
| `hwa` | float | First barrier curvature [MeV] |
| `Symb` | str | Second barrier symmetry designation |
| `Vb` | float | Second barrier height [MeV] |
| `hwb` | float | Second barrier curvature [MeV] |
| `Vc` | float | Third barrier height [MeV] |
| `Deltaf` | float | Pairing correction [MeV] |

### empirical_barriers  (77 entries, key=Nuclide, entry=Entry)

| Field | Type | Description |
|-------|------|-------------|
| `n` | Nuclide | Target nucleus |
| `s` | str | Element symbol |
| `syma` | str | Inner saddle-point symmetry (S/GA/MA) |
| `Va` | float | Inner barrier height [MeV] |
| `dVa` | float | Uncertainty on Va [MeV] |
| `hwa` | float | Inner barrier curvature [MeV] |
| `dwa` | float | Uncertainty on hwa [MeV] |
| `symb` | str | Outer saddle-point symmetry (S/GA/MA) |
| `Vb` | float | Outer barrier height [MeV] |
| `dVb` | float | Uncertainty on Vb [MeV] |
| `hwb` | float | Outer barrier curvature [MeV] |
| `dwb` | float | Uncertainty on hwb [MeV] |

### empirical_barriers_new  (77 entries, key=Nuclide, entry=Entry)

| Field | Type | Description |
|-------|------|-------------|
| `n` | Nuclide | Target nucleus |
| `s` | str | Element symbol |
| `syma` | str | Inner saddle-point symmetry (S/GA/MA) |
| `Va` | float | Inner barrier height [MeV] |
| `dVa` | float | Uncertainty on Va [MeV] |
| `hwa` | float | Inner barrier curvature [MeV] |
| `dwa` | float | Uncertainty on hwa [MeV] |
| `symb` | str | Outer saddle-point symmetry (S/GA/MA) |
| `Vb` | float | Outer barrier height [MeV] |
| `dVb` | float | Uncertainty on Vb [MeV] |
| `hwb` | float | Outer barrier curvature [MeV] |
| `dwb` | float | Uncertainty on hwb [MeV] |

### hfb_barriers  (0 entries, key=Nuclide, entry=None)

_(no field metadata)_

### hfbpath_bskg3
_legacy/absent in this release_

### hfbpath_d1m  (45 entries, key=Nuclide, entry=Entry)

| Field | Type | Description |
|-------|------|-------------|
| `n` | Nuclide | Target nucleus |
| `Z` | int | Charge number |
| `A` | int | Mass number |
| `sym` | str | Element symbol from the title line |
| `nbeta` | int | Number of beta deformation points |
| `Egs` | float | Ground-state binding energy [MeV] |
| `beta20` | list | Quadrupole deformation parameter |
| `beta22` | list | Quadrupole (triaxial) deformation parameter |
| `beta30` | list | Octupole deformation parameter |
| `dE` | list | Energy above ground state E-Egs [MeV] |
| `MI` | list | Collective inertia MI_ATD [hbar^2/MeV] (wrt line number) |
| `idx` | list | Index flag: -1 ground state, 1 barrier, 2 well |

### nld_fis  (2467 entries, key=Nuclide, entry=Entry)

| Field | Type | Description |
|-------|------|-------------|
| `n` | Nuclide | Target nucleus |
| `Z` | int | Charge number |
| `A` | int | Mass number |
| `label` | str | Saddle/well label (Max1/Max2/Max3/Min1/Min2) |
| `positive_parity` | riplpy.fission.nld_fis.ParityData | Positive-parity level density (ParityData) |
| `negative_parity` | riplpy.fission.nld_fis.ParityData | Negative-parity level density (ParityData) |

### rmf_axial  (45 entries, key=Nuclide, entry=Entry)

| Field | Type | Description |
|-------|------|-------------|
| `n` | Nuclide | Target nucleus |
| `Z` | int | Charge number |
| `A` | int | Mass number |
| `sym` | str | Element symbol from the header |
| `nbeta` | int | Number of beta deformation points |
| `Egs` | float | Ground-state energy [MeV] |
| `E0` | float | E(0+) band-head energy [MeV] |
| `beta20` | list | Quadrupole deformation parameter |
| `beta22` | list | Quadrupole (triaxial) deformation parameter |
| `beta30` | list | Octupole deformation parameter |
| `dE` | list | Energy relative to ground state E-EGS [MeV] |
| `Mu` | list | Collective inertia Mu |

### rmf_triaxial
_legacy/absent in this release_

## gamma

### experiment_mlo  (157 entries, key=Nuclide, entry=Entry)

| Field | Type | Description |
|-------|------|-------------|
| `n` | Nuclide | Target nucleus |
| `E1` | float | First GDR peak energy [MeV] |
| `W1` | float | First GDR peak width [MeV] |
| `E2` | float | Second GDR peak energy [MeV] |
| `W2` | float | Second GDR peak width [MeV] |
| `Er1` | float | First component energy [MeV] |
| `dEr1` | float | Uncertainty on Er1 [MeV] |
| `Wr1` | float | First component width [MeV] |
| `dWr1` | float | Uncertainty on Wr1 [MeV] |
| `Sr1` | float | First component strength [TRK units] |
| `dSr1` | float | Uncertainty on Sr1 |
| `CSp1` | float | First component Lorentzian peak cross section [mb] |
| `dCSp1` | float | Uncertainty on CSp1 [mb] |
| `Er2` | float | Second component energy [MeV] |
| `dEr2` | float | Uncertainty on Er2 [MeV] |
| `Wr2` | float | Second component width [MeV] |
| `dWr2` | float | Uncertainty on Wr2 [MeV] |
| `Sr2` | float | Second component strength [TRK units] |
| `dSr2` | float | Uncertainty on Sr2 |
| `CSp2` | float | Second component Lorentzian peak cross section [mb] |
| `dCSp2` | float | Uncertainty on CSp2 [mb] |
| `Sr` | float | Total strength Sr1+Sr2 [TRK units] |
| `dSr` | float | Uncertainty on Sr |
| `Id` | str | Identifier string |
| `Emin` | float | Lower energy limit of fit [MeV] |
| `Emax` | float | Upper energy limit of fit [MeV] |
| `reac` | int | Reaction flag |
| `reference` | str | Reference key |

### experiment_slo  (157 entries, key=Nuclide, entry=Entry)

| Field | Type | Description |
|-------|------|-------------|
| `n` | Nuclide | Target nucleus |
| `E1` | float | First GDR peak energy [MeV] |
| `W1` | float | First GDR peak width [MeV] |
| `E2` | float | Second GDR peak energy [MeV] |
| `W2` | float | Second GDR peak width [MeV] |
| `Er1` | float | First component energy [MeV] |
| `dEr1` | float | Uncertainty on Er1 [MeV] |
| `Wr1` | float | First component width [MeV] |
| `dWr1` | float | Uncertainty on Wr1 [MeV] |
| `Sr1` | float | First component strength [TRK units] |
| `dSr1` | float | Uncertainty on Sr1 |
| `CSp1` | float | First component Lorentzian peak cross section [mb] |
| `dCSp1` | float | Uncertainty on CSp1 [mb] |
| `Er2` | float | Second component energy [MeV] |
| `dEr2` | float | Uncertainty on Er2 [MeV] |
| `Wr2` | float | Second component width [MeV] |
| `dWr2` | float | Uncertainty on Wr2 [MeV] |
| `Sr2` | float | Second component strength [TRK units] |
| `dSr2` | float | Uncertainty on Sr2 |
| `CSp2` | float | Second component Lorentzian peak cross section [mb] |
| `dCSp2` | float | Uncertainty on CSp2 [mb] |
| `Sr` | float | Total strength Sr1+Sr2 [TRK units] |
| `dSr` | float | Uncertainty on Sr |
| `Id` | str | Identifier string |
| `Emin` | float | Lower energy limit of fit [MeV] |
| `Emax` | float | Upper energy limit of fit [MeV] |
| `reac` | int | Reaction flag |
| `reference` | str | Reference key |

### experiment_slo_errors  (157 entries, key=Nuclide, entry=Entry)

| Field | Type | Description |
|-------|------|-------------|
| `n` | Nuclide | Target nucleus |
| `E1` | float | First GDR peak energy [MeV] |
| `W1` | float | First GDR peak width [MeV] |
| `E2` | float | Second GDR peak energy [MeV] |
| `W2` | float | Second GDR peak width [MeV] |
| `Er1` | float | First component energy [MeV] |
| `dEr1` | float | Uncertainty on Er1 [MeV] |
| `Wr1` | float | First component width [MeV] |
| `dWr1` | float | Uncertainty on Wr1 [MeV] |
| `Sr1` | float | First component strength [TRK units] |
| `dSr1` | float | Uncertainty on Sr1 |
| `CSp1` | float | First component Lorentzian peak cross section [mb] |
| `dCSp1` | float | Uncertainty on CSp1 [mb] |
| `Er2` | float | Second component energy [MeV] |
| `dEr2` | float | Uncertainty on Er2 [MeV] |
| `Wr2` | float | Second component width [MeV] |
| `dWr2` | float | Uncertainty on Wr2 [MeV] |
| `Sr2` | float | Second component strength [TRK units] |
| `dSr2` | float | Uncertainty on Sr2 |
| `CSp2` | float | Second component Lorentzian peak cross section [mb] |
| `dCSp2` | float | Uncertainty on CSp2 [mb] |
| `Sr` | float | Total strength Sr1+Sr2 [TRK units] |
| `dSr` | float | Uncertainty on Sr |
| `Id` | str | Identifier string |
| `Emin` | float | Lower energy limit of fit [MeV] |
| `Emax` | float | Upper energy limit of fit [MeV] |
| `reac` | int | Reaction flag |
| `reference` | str | Reference key |

### experiment_smlo  (157 entries, key=Nuclide, entry=Entry)

| Field | Type | Description |
|-------|------|-------------|
| `n` | Nuclide | Target nucleus |
| `E1` | float | First GDR peak energy [MeV] |
| `W1` | float | First GDR peak width [MeV] |
| `E2` | float | Second GDR peak energy [MeV] |
| `W2` | float | Second GDR peak width [MeV] |
| `Er1` | float | First component energy [MeV] |
| `dEr1` | float | Uncertainty on Er1 [MeV] |
| `Wr1` | float | First component width [MeV] |
| `dWr1` | float | Uncertainty on Wr1 [MeV] |
| `Sr1` | float | First component strength [TRK units] |
| `dSr1` | float | Uncertainty on Sr1 |
| `CSp1` | float | First component Lorentzian peak cross section [mb] |
| `dCSp1` | float | Uncertainty on CSp1 [mb] |
| `Er2` | float | Second component energy [MeV] |
| `dEr2` | float | Uncertainty on Er2 [MeV] |
| `Wr2` | float | Second component width [MeV] |
| `dWr2` | float | Uncertainty on Wr2 [MeV] |
| `Sr2` | float | Second component strength [TRK units] |
| `dSr2` | float | Uncertainty on Sr2 |
| `CSp2` | float | Second component Lorentzian peak cross section [mb] |
| `dCSp2` | float | Uncertainty on CSp2 [mb] |
| `Sr` | float | Total strength Sr1+Sr2 [TRK units] |
| `dSr` | float | Uncertainty on Sr |
| `Id` | str | Identifier string |
| `Emin` | float | Lower energy limit of fit [MeV] |
| `Emax` | float | Upper energy limit of fit [MeV] |
| `reac` | int | Reaction flag |
| `reference` | str | Reference key |

### experiment_smlo_errors  (157 entries, key=Nuclide, entry=Entry)

| Field | Type | Description |
|-------|------|-------------|
| `n` | Nuclide | Target nucleus |
| `E1` | float | First GDR peak energy [MeV] |
| `W1` | float | First GDR peak width [MeV] |
| `E2` | float | Second GDR peak energy [MeV] |
| `W2` | float | Second GDR peak width [MeV] |
| `Er1` | float | First component energy [MeV] |
| `dEr1` | float | Uncertainty on Er1 [MeV] |
| `Wr1` | float | First component width [MeV] |
| `dWr1` | float | Uncertainty on Wr1 [MeV] |
| `Sr1` | float | First component strength [TRK units] |
| `dSr1` | float | Uncertainty on Sr1 |
| `CSp1` | float | First component Lorentzian peak cross section [mb] |
| `dCSp1` | float | Uncertainty on CSp1 [mb] |
| `Er2` | float | Second component energy [MeV] |
| `dEr2` | float | Uncertainty on Er2 [MeV] |
| `Wr2` | float | Second component width [MeV] |
| `dWr2` | float | Uncertainty on Wr2 [MeV] |
| `Sr2` | float | Second component strength [TRK units] |
| `dSr2` | float | Uncertainty on Sr2 |
| `CSp2` | float | Second component Lorentzian peak cross section [mb] |
| `dCSp2` | float | Uncertainty on CSp2 [mb] |
| `Sr` | float | Total strength Sr1+Sr2 [TRK units] |
| `dSr` | float | Uncertainty on Sr |
| `Id` | str | Identifier string |
| `Emin` | float | Lower energy limit of fit [MeV] |
| `Emax` | float | Upper energy limit of fit [MeV] |
| `reac` | int | Reaction flag |
| `reference` | str | Reference key |

### experiment_systematics_slo  (8980 entries, key=Nuclide, entry=Entry)

| Field | Type | Description |
|-------|------|-------------|
| `n` | Nuclide | Target nucleus |
| `E1` | float | First GDR peak energy [MeV] |
| `W1` | float | First GDR peak width [MeV] |
| `E2` | float | Second GDR peak energy [MeV] |
| `W2` | float | Second GDR peak width [MeV] |
| `Sr1` | float | First component strength [TRK units] |
| `CSp1` | float | First component Lorentzian peak cross section [mb] |
| `Sr2` | float | Second component strength [TRK units] |
| `CSp2` | float | Second component Lorentzian peak cross section [mb] |
| `Sr` | float | Total strength Sr1+Sr2 [TRK units] |
| `In` | int | Source flag (1 = experimental, 0 = systematics) |

### experiment_systematics_smlo  (8980 entries, key=Nuclide, entry=Entry)

| Field | Type | Description |
|-------|------|-------------|
| `n` | Nuclide | Target nucleus |
| `E1` | float | First GDR peak energy [MeV] |
| `W1` | float | First GDR peak width [MeV] |
| `E2` | float | Second GDR peak energy [MeV] |
| `W2` | float | Second GDR peak width [MeV] |
| `Sr1` | float | First component strength [TRK units] |
| `CSp1` | float | First component Lorentzian peak cross section [mb] |
| `Sr2` | float | Second component strength [TRK units] |
| `CSp2` | float | Second component Lorentzian peak cross section [mb] |
| `Sr` | float | Total strength Sr1+Sr2 [TRK units] |
| `In` | int | Source flag (1 = experimental, 0 = systematics) |

### gsf  (1 entries, key=Nuclide, entry=PacketEntry)

| Field | Type | Description |
|-------|------|-------------|
| `U` | list | Photon-energy grid [MeV] |
| `fE1` | list | E1 strength column [mb/MeV for legacy gamma-strength-micro files; MeV^-3 for the D1M+QRPA fallback] |

### gsf_d1m  (43 entries, key=Nuclide, entry=Entry)

| Field | Type | Description |
|-------|------|-------------|
| `n` | Nuclide | Target nucleus |
| `U` | list | Photon-energy grid [MeV] |
| `fE1` | list | Cold-nucleus (T=0) E1 photon strength function [MeV^-3] |
| `T` | list | Temperature/excitation column labels [MeV] |
| `fE1_T` | list | Per-temperature E1 strength rows (rows = U, cols = T) [MeV^-3] |

### psf  (103 entries, key=Nuclide, entry=Entry)

| Field | Type | Description |
|-------|------|-------------|
| `n` | Nuclide | Target nucleus |
| `E` | list | Photon energy grid [MeV] |
| `f` | list | Photon strength function [MeV^-3] |
| `rows` | list | Raw numeric data table (one list per data line) |
| `columns` | list | Column labels parsed from the file header |
| `source` | str | Experimental technique category (subdirectory name) |
| `filename` | str | Source data file name |

### smlo_e1  (1 entries, key=Nuclide, entry=Entry)

| Field | Type | Description |
|-------|------|-------------|
| `n` | Nuclide | Target nucleus |
| `U` | list | Photon-energy grid [MeV] |
| `fE1` | list | Cold-nucleus (T=0) E1 photon strength function [MeV^-3] |
| `fE1_T` | list | Per-temperature E1 strength rows (rows = U, cols = T) [MeV^-3] |
| `T` | list | Temperature column labels (T=0.0..2.0 MeV) |
| `Er1` | float | SMLO Lorentzian peak energy [MeV] |
| `Wr1` | float | SMLO Lorentzian width [MeV] |
| `S1` | float | SMLO Lorentzian strength [TRK units] |
| `beta` | float | Quadrupole deformation parameter |

### smlo_m1  (43 entries, key=Nuclide, entry=Entry)

| Field | Type | Description |
|-------|------|-------------|
| `n` | Nuclide | Target nucleus |
| `U` | list | Photon-energy grid [MeV] |
| `fE1` | list | Cold-nucleus (T=0) M1 photon strength function [MeV^-3] |
| `T` | list | Strength column labels (T=0, T>0) |
| `fE1_T` | list | M1 strength rows (rows = U, cols = T=0 / T>0) [MeV^-3] |

### theory_gdr  (0 entries, key=Nuclide, entry=None)

_(no field metadata)_

### tlo  (53 entries, key=Nuclide, entry=Entry)

| Field | Type | Description |
|-------|------|-------------|
| `n` | Nuclide | Target nucleus |
| `U` | list | Photon-energy grid [MeV] |
| `fE1` | list | E1 dipole strength function [MeV^-3] |
| `fE1_T` | list | Full per-nucleus strength rows (rows = U) [MeV^-3] |
| `beta` | float | Ground-state beta deformation parameter (Delaroche et al.) |
| `gamma` | float | Ground-state gamma deformation parameter [deg] |
| `mode` | str | Deformation mode flag (EFF = reduced near shells, INTER = interpolated odd nuclei, 5DCH = Delaroche et al.) |

## levels

### constant_temperature  (3557 entries, key=Nuclide, entry=Entry)

| Field | Type | Description |
|-------|------|-------------|
| `n` | Nuclide | Target nucleus |
| `T` | float | Nuclear temperature [MeV] |
| `dT` | float | Uncertainty on T [MeV] |
| `U0` | float | Back-shift energy [MeV] |
| `dU0` | float | Uncertainty on U0 [MeV] |
| `Nlev` | int | Total number of ENSDF levels |
| `Nmax` | int | Completeness limit (level number) |
| `N0` | int | First level included in fit |
| `Nc` | int | Level where unique spin sequence ends |
| `Umax` | float | Energy at completeness limit [MeV] |
| `Uc` | float | Energy at Nc [MeV] |
| `Chi` | float | Fit quality chi-squared |
| `Fit` | str | Fit flag (blank if no fit) |
| `Flag` | str | F if poor fit (Chi>0.05) |
| `NoX` | int | Number of levels with +X notation |
| `Xm` | int | First level with +X notation |
| `Ex` | float | Energy of first +X level [MeV] |
| `sigma` | float | Spin cutoff parameter |

### discrete_levels  (3557 entries, key=Nuclide, entry=LevelEntry)

| Field | Type | Description |
|-------|------|-------------|
| `nucleus` | Nuclide | Target nucleus |
| `levels` | list | List of discrete levels (ordered by energy) |
| `gammas` | dict | Dict of gamma transitions by level index |
| `extras` | list | Additional ENSDF metadata |

## masses

### ame20  (3558 entries, key=Nuclide, entry=Entry)

| Field | Type | Description |
|-------|------|-------------|
| `n` | Nuclide | Target nucleus |
| `flag` | int | Data quality flag (0=measured, 1=extrapolated) |
| `Mexp` | float | Experimental mass excess [MeV] |
| `Err` | float | Uncertainty on mass excess [MeV] |

### bskg3  (8587 entries, key=Nuclide, entry=Entry)

| Field | Type | Description |
|-------|------|-------------|
| `n` | Nuclide | Target nucleus |
| `flag` | int | Data quality flag |
| `Mexp` | float | Experimental mass excess [MeV] |
| `Err` | float | Uncertainty on mass excess [MeV] |
| `Mth` | float | Calculated BSkG3 mass excess [MeV] |
| `beta20` | float | Quadrupole deformation (axial) |
| `beta22` | float | Quadrupole deformation (triaxial) |
| `Gamma` | float | Triaxial deformation angle [deg] |
| `beta30` | float | Octupole deformation (axial) |
| `beta32` | float | Octupole deformation (triaxial) |
| `beta4` | float | Hexadecapole deformation |
| `Rch` | float | Charge radius [fm] |

### d1m  (8447 entries, key=Nuclide, entry=Entry)

| Field | Type | Description |
|-------|------|-------------|
| `n` | Nuclide | Target nucleus |
| `flag` | int | Data quality flag |
| `Mexp` | float | Experimental mass excess [MeV] |
| `Err` | float | Uncertainty on mass excess [MeV] |
| `Mth` | float | Calculated D1M mass excess [MeV] |
| `beta20` | float | Quadrupole deformation |
| `Rch` | float | Charge radius [fm] |

### deformations  (328 entries, key=Nuclide, entry=Entry)

| Field | Type | Description |
|-------|------|-------------|
| `n` | Nuclide | Target nucleus |
| `beta2` | float | Experimental quadrupole deformation |
| `error` | float | Uncertainty on beta2 |

### density_bskg3
_legacy/absent in this release_

### density_d1m
_legacy/absent in this release_

### frdm1995  (0 entries, key=Nuclide, entry=None)

_(no field metadata)_

### frdm2012  (9420 entries, key=Nuclide, entry=Entry)

| Field | Type | Description |
|-------|------|-------------|
| `n` | Nuclide | Target nucleus |
| `flag` | int | Data quality flag (0=no exp, 1=recommended, 2=measured) |
| `Mexp` | float | Experimental mass excess [MeV] |
| `Err` | float | Uncertainty on mass excess [MeV] |
| `Mth` | float | Calculated FRDM2012 mass excess [MeV] |
| `Emic` | float | Microscopic correction energy [MeV] |
| `beta2` | float | Quadrupole deformation |
| `beta3` | float | Octupole deformation |
| `beta4` | float | Hexadecapole deformation |
| `beta6` | float | Hexacontatetrapole deformation |

### hfb14  (0 entries, key=Nuclide, entry=None)

_(no field metadata)_

### hfb27  (10277 entries, key=Nuclide, entry=Entry)

| Field | Type | Description |
|-------|------|-------------|
| `n` | Nuclide | Target nucleus |
| `flag` | int | Data quality flag |
| `Mexp` | float | Experimental mass excess [MeV] |
| `Err` | float | Uncertainty on mass excess [MeV] |
| `Mth` | float | Calculated HFB-27 mass excess [MeV] |
| `beta20` | float | Quadrupole deformation |
| `beta40` | float | Hexadecapole deformation |
| `Rch` | float | Charge radius [fm] |

### natab  (289 entries, key=Nuclide, entry=Entry)

| Field | Type | Description |
|-------|------|-------------|
| `n` | Nuclide | Target nucleus |
| `abundance` | float | Natural isotopic abundance [%] |
| `uncertainty` | float | Uncertainty on abundance [%] |

### ws4  (10339 entries, key=Nuclide, entry=Entry)

| Field | Type | Description |
|-------|------|-------------|
| `n` | Nuclide | Target nucleus |
| `flag` | int | Data quality flag |
| `Mexp` | float | Experimental mass excess [MeV] |
| `Err` | float | Uncertainty on mass excess [MeV] |
| `Mth` | float | Calculated WS4 mass excess [MeV] |
| `Esh` | float | Shell correction energy [MeV] |
| `beta2` | float | Quadrupole deformation |
| `beta4` | float | Hexadecapole deformation |
| `beta6` | float | Hexacontatetrapole deformation |

## optical

### deformations  (1423 entries, key=tuple, entry=DeformationEntry)

_(no field metadata)_

### index  (581 entries, key=int, entry=IndexEntry)

_(no field metadata)_

### potentials  (584 entries, key=int, entry=CoupledChannelOMP)

_(no field metadata)_

### references  (143 entries, key=int, entry=Reference)

_(no field metadata)_

### rop2013  (79 entries, key=tuple, entry=ROP2013Entry)

_(no field metadata)_

## resonances

### pwave  (248 entries, key=Nuclide, entry=Entry)

| Field | Type | Description |
|-------|------|-------------|
| `n` | Nuclide | Target nucleus |
| `sym` | str | Element symbol of target nucleus |
| `L` | int | Angular momentum of incident neutron (0=s-wave, 1=p-wave) |
| `Io` | float | Spin of the target ground state |
| `parity` | str | Parity of the target ground state |
| `Bn` | float | Neutron binding energy [MeV] |
| `D` | float | Average resonance spacing (RIPL-3) [eV] |
| `Derr` | float | Uncertainty on D (RIPL-3) [eV] |
| `Gam` | float | Average radiative width (RIPL-3) [eV] |
| `Gerr` | float | Uncertainty on Gam (RIPL-3) [eV] |
| `Str` | float | Neutron strength function (RIPL-3) [10^-4] |
| `Serr` | float | Uncertainty on Str (RIPL-3) [10^-4] |
| `D_BNL` | float | Average resonance spacing (BNL) [eV] |
| `D_BNL_err` | float | Uncertainty on D_BNL [eV] |
| `Gam_BNL` | float | Average radiative width (BNL) [eV] |
| `Gam_BNL_err` | float | Uncertainty on Gam_BNL [eV] |
| `Str_BNL` | float | Neutron strength function (BNL) [10^-4] |
| `Str_BNL_err` | float | Uncertainty on Str_BNL [10^-4] |

### swave  (324 entries, key=Nuclide, entry=Entry)

| Field | Type | Description |
|-------|------|-------------|
| `n` | Nuclide | Target nucleus |
| `sym` | str | Element symbol of target nucleus |
| `L` | int | Angular momentum of incident neutron (0=s-wave, 1=p-wave) |
| `Io` | float | Spin of the target ground state |
| `parity` | str | Parity of the target ground state |
| `Bn` | float | Neutron binding energy [MeV] |
| `D` | float | Average resonance spacing (RIPL-3) [eV] |
| `Derr` | float | Uncertainty on D (RIPL-3) [eV] |
| `Gam` | float | Average radiative width (RIPL-3) [eV] |
| `Gerr` | float | Uncertainty on Gam (RIPL-3) [eV] |
| `Str` | float | Neutron strength function (RIPL-3) [10^-4] |
| `Serr` | float | Uncertainty on Str (RIPL-3) [10^-4] |
| `D_BNL` | float | Average resonance spacing (BNL) [eV] |
| `D_BNL_err` | float | Uncertainty on D_BNL [eV] |
| `Gam_BNL` | float | Average radiative width (BNL) [eV] |
| `Gam_BNL_err` | float | Uncertainty on Gam_BNL [eV] |
| `Str_BNL` | float | Neutron strength function (BNL) [10^-4] |
| `Str_BNL_err` | float | Uncertainty on Str_BNL [10^-4] |
