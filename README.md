# RIPL-4

The latest RIPL Coordinated Research Project (CRP) on "Parameters for Calculation of Nuclear Reactions of Relevance to Non-Energy Nuclear Applications" - RIPL-4 was successfully completed in 2024, after many years of challenging work carried out through four consecutive IAEA coordinated research projects.
The RIPL-4 library was released in May 2026, and is available at https://nds.iaea.org/RIPL-4/ as well as in this IAEA GitHub repository. The numerical data and computer codes included in RIPL-4 are arranged in seven segments (or directories):

MASSES: contains the experimental and recommended masses evaluated by Wang et al. (2021) as well as the ground-state properties for about 9000 nuclei obtained within five different mass models, two of which are based on the macroscopic-microscopic approach and three on self-consistent mean-field models.

LEVELS: contains 118 datasets (one for each element) with all known level schemes and electromagnetic decay probabilities available from ENSDF as of October 2025 complemented by information from the NUBASE 2020 database.

RESONANCES: contains average resonance parameters obtained from both RIPL-3 and Brookhaven National Laboratory for 324 nuclei.

OPTICAL: contains 495 sets of phenomenological optical model parameters defined in a wide energy range. When there is insufficient experimental data to define the parameters of the optical phenomenological model, the evaluator has to resort to either global parameterizations or microscopic
approaches. Radial density distributions to be used as input for microscopic optical-model calculations are stored in the MASSES segment.

DENSITIES: contains phenomenological parameterizations based on the modified Fermi gas model similar to the RIPL-3 recommendations, but also updated microscopic calculations based on the combinatorial approach with realistic microscopic single-particle level schemes. All recommended level density models are consistent with the average neutron resonance spacings, discrete level schemes, and Oslo data.

GAMMA: contains GDR resonance parameters as well as large-scale photon strength function calculations based on the phenomenological Lorentzian approach and microscopic mean-field (D1M) plus QRPA recommended by the CRP on Updating the Photonuclear Data Library and Generating a Reference Database for Photon Strength Functions (2016--2019). The data can also be accessed from https://nds.iaea.org/PSFdatabase. Published in S. Goriely et al., The European Physical Journal A 55, 172 (2019), 10.1140/epja/i2019-12840-1. The GDR parameters were published in V.A. Plujko, et al., At. Data Nucl. Data Tables 123, 1 (2018), 10.1016/j.adt.2018.03.002.

FISSION: includes prescriptions for optimized calculations of the neutron-induced fission cross section for known actinides, including in particular a revised set of empirical barriers. This segment also features global prescriptions for fission paths with associated collective inertias and nuclear level densities at fission saddle points based from different sources: macroscopic-microscopic approaches and both non-relativistic (Skyrme and Gogny-type) and relativistic mean-field models.
