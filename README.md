# ignite

ignite suggests how to make an inorganic compound by solid-state synthesis. Type a formula at https://material-hub.github.io/ignite/ and it returns precursor sets, a firing temperature with a range, the furnace atmosphere and the impurity phases seen in similar syntheses.

This demo copy holds the stored plans of 20 compounds. The full library of 13,031 compounds is on the live site. For each compound in the library the full plan was computed ahead of time. It adds Materials Project reaction energies, NIST-JANAF heating data, crucible checks, oxygen vacancies, ReactCA phase evolution and phonons. Any other formula gets a lighter plan, computed in your browser from 64,375 published recipes of 30,993 compounds. That plan uses the same learned models and leaves out the physics.

## Data collected with language models

ignite is built on synthesis data that text mining and language models read out of published papers.

- The 64,375 recipes come from the text-mined dataset of Kononova et al. (2019) and the dataset of Lee et al. (2025), which was mined with a large language model. Lee et al. also give the impurity labels.
- Our own language-model agents, mostly Claude, read 1,233 newer papers and extracted 3,008 synthesis records across 12 materials families. Each record was checked against its paper: 1,135 passed, 1,739 wait for expert review and 134 were rejected. These records are kept out of the planner's knowledge base. The Claude benchmark questions come from the papers published in 2025 or later.
- Claude writes the optional lab protocol, and the Valency MCP was used to check related work.

## How it works

Precursor sets come from the recipes of similar known compounds and are ranked by a model trained on the Kononova et al. recipes. The firing temperature is a Huber regression on the nearest recipes and on element properties, with a conformal range. When nothing in the data is close enough, it gives no temperature and says so. The atmosphere comes from a small classifier. The dwell time and heating steps come from similar recipes. Impurity phases come from the Lee et al. 2025 labels of the nearest syntheses. Fixed rules mark each claim trust, caution, abstain or unavailable. None of this uses a language model. The Claude protocol tab is optional and calls Claude from your browser with your own API key.

## Accuracy on held-out data

| Task | Test items | ignite | Simple baseline |
|---|---|---|---|
| Right precursor set in top 1 | 805 targets | 65.1% | 59.6% (frequency prior) |
| Right precursor set in top 5 | 805 targets | 81.6% | 75.2% (frequency prior) |
| Right precursor set in top 5, A-Lab | 42 targets | 71.4% | 61.9% (frequency prior) |
| Firing temperature, mean error | 664 syntheses | 153.6 °C | 158.6 °C (neighbour median) |
| 90% range holds the true temperature | 664 syntheses | 92.9% | |
| Furnace atmosphere | 691 recipes | 80.5% | 76.6% (always air) |
| Reported impurity in top 3 | 958 queries | 51.0% | 32.9% (global frequency) |

The temperature model and the neighbour median are close. The 95% interval of their gap runs from -0.1 to 10.1 °C.

## Claude benchmark

In October 2026 we asked Claude (claude-opus-5-5) the same questions blind, on recipes published in 2025 or later. In the "Claude with ignite" column each question also carried ignite's own answer and its held-out accuracy.

| Task | Items | Claude | Claude with ignite | ignite |
|---|---|---|---|---|
| Right precursor set in top 5 | 242 | 53.7% | 71.5% | 71.5% |
| Right precursor set in top 5, papers after Claude's cutoff | 23 | 60.9% | 69.6% | 69.6% |
| Solid-state or sol-gel route | 240 | 82.1% | 82.1% | 62.9% |
| Furnace atmosphere | 159 | 95.0% | 96.2% | 82.4% |
| Dwell time, mean error | 261 | 4.2 h | 4.2 h | 4.7 h |
| Right precursor set in top 5, older text-mined set | 805 | 81.7% | 82.9% | 81.6% |

Claude's training data runs to June 2026, so most of these papers could be in it. Only the 23 papers published after that date are unseen, and on those every gap falls inside its 95% interval. With ignite's answer in front of it, Claude beat ignite alone on the route (82.1% against 62.9%), the atmosphere (96.2% against 82.4%) and the dwell time (4.2 h against 4.7 h), and each of these gaps lies outside its 95% interval. This table keeps the tasks where Claude with ignite scored as well as or better than both. On firing temperature Claude alone had the lower error (105.4 °C against 124.8 °C for ignite on 241 syntheses).

## Files

`index.html` is the page. `ignite_data/` holds 20 stored plans: SrTiO3, LiCoO2, CdSe0.5Te0.5, BaTiO3, BaZrS3, CaTiO3, CdTe, LiFePO4, LiNiO2, MgAl2O4, SrFeO3, BaZrO3, NaNbO3, KNbO3, LaCoO3, LaFeO3, LaAlO3, SrZrO3, CaZrO3, BaSnO3. `ignite_lite.js` is the in-browser planner. Its knowledge file is left out of this copy, so a formula outside the 20 is reported as not in the library. `version.json` lets the page reload itself after an update.

## Related work checked with Valency

On 22 October 2026 we searched the Valency literature index (arXiv, PubMed, ChemRxiv and other preprint servers) with three questions.

1. Machine learning that recommends precursors and conditions for solid-state synthesis, work from 2024 onward.
2. Language-model agents that plan inorganic synthesis from literature evidence, work from 2024 onward.
3. Conformal ranges or abstention for a predicted synthesis temperature, work from 2023 onward.

The first two searches returned 20 distinct works. 18 of them were already in our related-work survey of 6 October 2026, which was built from OpenAlex and arXiv. Two were new to us:

- Inada, Y. et al. Elemental Reactivity Maps for Materials Discovery. Chemistry of Materials 37, 2097-2105 (2025). https://doi.org/10.1021/acs.chemmater.4c02259
- Takahara, I., Mizoguchi, T. and Liu, B. Accelerated inorganic materials design with generative AI agents. Cell Reports Physical Science 6, 103019 (2025). https://doi.org/10.1016/j.xcrp.2025.103019

The third search returned conformal methods for other properties (steel fatigue strength, void swelling, polymer glass transition, interatomic potentials). None of the returned papers puts a range on a predicted firing temperature, which is what the ignite temperature model reports.

## References

### Recipe and outcome datasets

1. Kononova, O. et al. Text-mined dataset of inorganic materials synthesis recipes. Scientific Data 6, 203 (2019). https://doi.org/10.1038/s41597-019-0224-1. Base solid-state recipes, the ranker labels and the fixed split.
2. Lee, S. et al. Text-mined dataset of solid-state syntheses with impurity phases using Large Language Model. Scientific Data 12, 1969 (2025). https://doi.org/10.1038/s41597-025-06222-y. Added retrieval recipes and the impurity labels.
3. Chen, W. Y. et al. Text-mined dataset of sol-gel synthesis protocols with phase purity outcomes using large language models. ChemRxiv preprint (2026). https://doi.org/10.26434/chemrxiv.15006472/v1. Kept with the data for reference. It is a preprint and plays no part in the plans.
4. Szymanski, N. J. et al. An autonomous laboratory for the accelerated synthesis of inorganic materials. Nature 624, 86-91 (2023). https://doi.org/10.1038/s41586-023-06734-w. A-Lab robotic outcomes and the reacted or inert precursor pairs of its supplement.
5. Szymanski, N. J. et al. Author Correction: An autonomous laboratory for the accelerated synthesis of inorganic materials. Nature 650, E1 (2026). https://doi.org/10.1038/s41586-025-09992-y. The correction applied to the A-Lab outcome labels.
6. Szymanski, N. J. et al. Autonomous and dynamic precursor selection for solid-state materials synthesis. Nature Communications 14, 6956 (2023). https://doi.org/10.1038/s41467-023-42329-9. ARROWS outcomes. Only the experimentally verified runs are used.
7. Walters, L. et al. Precursor Genome. Zenodo, version 1.0 (2026). https://doi.org/10.5281/zenodo.21285546. Heated precursor pairs, the labels of the pairwise reaction model.

### Thermochemistry and oxygen vacancies

8. Jain, A. et al. Commentary: The Materials Project: A materials genome approach to accelerating materials innovation. APL Materials 1, 011002 (2013). https://doi.org/10.1063/1.4812323. GGA/GGA+U formation energies and hull distances for the reaction energies and the scorecard.
9. Chase, M. W., Jr. NIST-JANAF Thermochemical Tables, 4th edn. American Chemical Society and American Institute of Physics for the National Institute of Standards and Technology, Washington, DC and New York (1998). Measured Gibbs energies for thermal onsets and gas equilibria.
10. Kumagai, Y. et al. Insights into oxygen vacancies from high-throughput first-principles calculations. Physical Review Materials 5, 123803 (2021). https://doi.org/10.1103/PhysRevMaterials.5.123803. Oxygen-vacancy energies, charge states and band gaps.
11. Witman, M. D. et al. Defect graph neural networks for materials discovery in high-temperature clean-energy applications. Nature Computational Science 3, 675-686 (2023). https://doi.org/10.1038/s43588-023-00495-2. A second vacancy set, used only for the uncertainty of the vacancy energies.
12. Bale, C. W. et al. FactSage thermochemical software and databases. Calphad 26, 189-228 (2002). https://doi.org/10.1016/S0364-5916(02)00035-4. Gibbs energy of CO<sub>2</sub> from the FactSage table shipped with the reaction network.
13. Sun, W. et al. The thermodynamic scale of inorganic crystalline metastability. Science Advances 2, e1600225 (2016). https://doi.org/10.1126/sciadv.1600225. The metastability scale behind the energy-above-hull check.
14. Aykol, M. et al. Thermodynamic limit for synthesis of metastable inorganic materials. Science Advances 4, eaaq0148 (2018). https://doi.org/10.1126/sciadv.aaq0148. The amorphous limit behind the same check.

### Reference sets for compound queries

15. Gharakhanyan, V. et al. Discovering melting temperature prediction models of inorganic solids by combining supervised and unsupervised learning. The Journal of Chemical Physics 160, 204112 (2024). https://doi.org/10.1063/5.0207033. Melting points of binary compounds.
16. Lee, S. et al. Data-driven analysis of text-mined seed-mediated syntheses of gold nanoparticles. Digital Discovery 4, 93-104 (2025). https://doi.org/10.1039/d4dd00158c. Seed-mediated gold nanoparticle recipes.
17. Therrien, F. et al. OBELiX: a curated dataset of crystal structures and experimentally measured ionic conductivities for lithium solid-state electrolytes. Digital Discovery 5, 910-918 (2026). https://doi.org/10.1039/d5dd00441a. Measured lithium-ion conductivities.
18. Wang, Z. et al. Dataset of solution-based inorganic materials synthesis procedures extracted from the scientific literature. Scientific Data 9, 231 (2022). https://doi.org/10.1038/s41597-022-01317-2. Solution-based hydrothermal and precipitation routes.
19. Pan, E. et al. ZeoSyn: A Comprehensive Zeolite Synthesis Dataset Enabling Machine-Learning Rationalization of Hydrothermal Parameters. ACS Central Science 10, 729-743 (2024). https://doi.org/10.1021/acscentsci.3c01615. Zeolite syntheses.
20. Kim, S. et al. PubChem 2023 update. Nucleic Acids Research 51, D1373-D1380 (2023). https://doi.org/10.1093/nar/gkac956. Handbook values, precursor names and hazard statements.

### Synthesis studies used for comparison

21. Huo, H. et al. Machine-Learning Rationalization and Prediction of Solid-State Synthesis Conditions. Chemistry of Materials 34, 7323-7336 (2022). https://doi.org/10.1021/acs.chemmater.2c01293. Melting point and precursor stability as descriptors of the firing temperature.
22. Leeman, J. et al. Challenges in High-Throughput Inorganic Materials Prediction and Autonomous Synthesis. PRX Energy 3, 011002 (2024). https://doi.org/10.1103/PRXEnergy.3.011002. Limits of automated phase assignment in the A-Lab outcomes.
23. Tunç Parlak, T. et al. Formation of SrTiO<sub>3</sub> in mechanically activated SrCO<sub>3</sub>-TiO<sub>2</sub> system. Journal of Thermal Analysis and Calorimetry 127, 63-69 (2016). https://doi.org/10.1007/s10973-016-5385-y. SrTiO<sub>3</sub> from SrCO<sub>3</sub> and TiO<sub>2</sub> without mechanical activation.
24. Rocha-Rangel, E. et al. Synthesis of SrTiO<sub>3</sub> by the Calcination of SrCO<sub>3</sub> and TiO<sub>2</sub> Mixtures Intensively Ground by Means of High Energy Milling. Archives of Metallurgy and Materials 65, 621-626 (2020). https://doi.org/10.24425/amm.2020.132801. The SrCO<sub>3</sub> decomposition range and single-phase SrTiO<sub>3</sub> after high-energy milling.

### Simulation, models and software

25. Gallant, M. C. et al. A Cellular Automaton Simulation for Predicting Phase Evolution in Solid-State Reactions. Chemistry of Materials 37, 210-223 (2025). https://doi.org/10.1021/acs.chemmater.4c02301. The ReactCA cellular automaton behind the phase evolution tab.
26. McDermott, M. J. et al. A graph-based network for predicting chemical reaction pathways in solid-state materials synthesis. Nature Communications 12, 3097 (2021). https://doi.org/10.1038/s41467-021-23339-x. Enumeration of the reactions in each chemical system.
27. Bartel, C. J. et al. Physical descriptor for the Gibbs energy of inorganic crystalline solids and temperature-dependent materials chemistry. Nature Communications 9, 4168 (2018). https://doi.org/10.1038/s41467-018-06682-4. The Gibbs energy descriptor added to the formation enthalpies.
28. Batatia, I. et al. A foundation model for atomistic materials chemistry. The Journal of Chemical Physics 163, 184110 (2025). https://doi.org/10.1063/5.0297006. The MACE-MPA-0 model used for the harmonic phonons.
29. Togo, A. First-principles Phonon Calculations with Phonopy and Phono3py. Journal of the Physical Society of Japan 92, 012001 (2023). https://doi.org/10.7566/jpsj.92.012001. Phonopy phonon calculations.
30. Ong, S. P. et al. Python Materials Genomics (pymatgen): A robust, open-source python library for materials analysis. Computational Materials Science 68, 314-319 (2013). https://doi.org/10.1016/j.commatsci.2012.10.028. Formula parsing and oxidation-state guesses.

### Statistical methods

31. Byrd, R. H. et al. A Limited Memory Algorithm for Bound Constrained Optimization. SIAM Journal on Scientific Computing 16, 1190-1208 (1995). https://doi.org/10.1137/0916069. The L-BFGS-B fit of the precursor ranker.
32. Huber, P. J. Robust Estimation of a Location Parameter. The Annals of Mathematical Statistics 35, 73-101 (1964). https://doi.org/10.1214/aoms/1177703732. Huber regression in the firing-temperature model.
33. Efron, B. Bootstrap Methods: Another Look at the Jackknife. The Annals of Statistics 7, 1-26 (1979). https://doi.org/10.1214/aos/1176344552. Bootstrap intervals of the test metrics.
34. Hanley, J. A. & McNeil, B. J. The meaning and use of the area under a receiver operating characteristic (ROC) curve. Radiology 143, 29-36 (1982). https://doi.org/10.1148/radiology.143.1.7063747. The ROC area of the outcome tests.

### Energy materials used as examples

35. Alves dos Santos, L. et al. Comparison of Synthesized Microstructured and Commercial FePO<sub>4</sub> as Precursors for High-Performance LiFePO<sub>4</sub>/C Cathode Materials. ACS Omega 11, 43682-43693 (2026). https://doi.org/10.1021/acsomega.6c02734. LiFePO<sub>4</sub>, lithium-ion battery cathode.
36. Pandit, B. et al. Na<sub>3</sub>V<sub>2</sub>(PO<sub>4</sub>)<sub>3</sub> Cathode for Room-Temperature Solid-State Sodium-Ion Batteries: Advanced In Situ Synchrotron X-ray Studies to Understand Intermediate Phase Evolution. Chemistry of Materials 36, 2314-2324 (2024). https://doi.org/10.1021/acs.chemmater.3c02585. Na<sub>3</sub>V<sub>2</sub>(PO<sub>4</sub>)<sub>3</sub>, sodium-ion battery cathode.
37. Cheng, E. J. et al. Degradation Mechanism of All-Solid-State Li-Metal Batteries Studied by Electrochemical Impedance Spectroscopy. ACS Applied Materials & Interfaces 14, 40881-40889 (2022). https://doi.org/10.1021/acsami.2c09841. Li<sub>7</sub>La<sub>3</sub>Zr<sub>2</sub>O<sub>12</sub>, oxide solid electrolyte.
38. Toan, T. V. et al. AC conductivity study of mechanochemically synthesized solid electrolytes of Li<sub>6-a</sub>M<sup>n</sup><sub>a/n</sub>PS<sub>5</sub>Cl (M = Ca, Mg, Ba, Zn, Al, Y). RSC Advances 14, 16240-16247 (2024). https://doi.org/10.1039/d4ra02621g. Li<sub>6</sub>PS<sub>5</sub>Cl, sulfide solid electrolyte.
39. van der Maas, E. et al. Re-investigating the structure-property relationship of the solid electrolytes Li<sub>3-x</sub>In<sub>1-x</sub>Zr<sub>x</sub>Cl<sub>6</sub> and the impact of In-Zr(IV) substitution. Journal of Materials Chemistry A 11, 4559-4571 (2023). https://doi.org/10.1039/d2ta08433c. Li<sub>3</sub>InCl<sub>6</sub>, halide solid electrolyte.
40. Rodkey, N. et al. Pulsed Laser Deposition of Cs<sub>2</sub>AgBiBr<sub>6</sub>: from Mechanochemically Synthesized Powders to Dry, Single-Step Deposition. Chemistry of Materials 33, 7417-7422 (2021). https://doi.org/10.1021/acs.chemmater.1c02054. Cs<sub>2</sub>AgBiBr<sub>6</sub>, lead-free halide solar absorber.
41. Riva, S. et al. Electronic Structure and Surface Chemistry of BaZrS<sub>3</sub> Perovskite Powder and Sputtered Thin Film. ACS Applied Materials & Interfaces 16, 40210-40221 (2024). https://doi.org/10.1021/acsami.4c06758. BaZrS<sub>3</sub>, chalcogenide solar absorber.
42. Zhang, S. et al. Achieving enhanced and robust propane selective oxidative dehydrogenation via electrically-driven continuous chemical looping. Nature Communications 17, 9084 (2026). https://doi.org/10.1038/s41467-026-76063-9. La<sub>0.6</sub>Sr<sub>0.4</sub>Co<sub>0.2</sub>Fe<sub>0.8</sub>O<sub>3</sub>, oxygen-variable catalyst and fuel-cell oxide.
43. Luo, Y. et al. Spatiotemporal alignment of hole transfer and water oxidation for highly efficient photocatalytic water splitting. Nature Communications 17, 2767 (2026). https://doi.org/10.1038/s41467-026-69276-5. SrTiO<sub>3</sub>, oxide photocatalyst.
44. Jung, S.-H. et al. Microstructural Engineering for Power Factor and ZT Optimization in Cu Alloyed n-Type Bi<sub>2-x</sub>Sb<sub>x</sub>Te<sub>3</sub> Thermoelectric Materials. ACS Omega (2025). https://doi.org/10.1021/acsomega.5c07507. Bi<sub>1.7</sub>Sb<sub>0.3</sub>Te<sub>3</sub>, thermoelectric.
