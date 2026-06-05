# Teaching Module: Computational Spectroscopy Labs

Three independent graduate-level lab modules built on the Q-Chem workflow
in the parent repository. Sessions can be run in any order or combination;
each is roughly 3 hours and self-contained.

## Sessions

| File | Spectroscopy | Molecule | Focus |
|------|--------------|----------|-------|
| [`Session1_UVVis_Benchmarking_Handout.docx`](Session1_UVVis_Benchmarking_Handout.docx) | UV-Vis (TDDFT) | para-nitroaniline | Functional and basis-set sensitivity for charge-transfer states |
| [`Session2_ECD_Conformational_Averaging_Handout.docx`](Session2_ECD_Conformational_Averaging_Handout.docx) | ECD (TDDFT) | (1R,2S)-ephedrine | Conformer search, Boltzmann averaging, and absolute-configuration assignment |
| [`Session3_VCD_IR_Handout.docx`](Session3_VCD_IR_Handout.docx) | IR and VCD (analytic Hessian) | (R)-methyloxirane | Analytic frequency calculation, VCD sign patterns, and frequency scaling |

## Other materials

- [`Instructor_Lesson_Plan.docx`](Instructor_Lesson_Plan.docx): module
  overview, prerequisites, software requirements, suggested timing, common
  student difficulties, and assessment guidance.
- [`Student_Report_Template_UVVis_ECD.docx`](Student_Report_Template_UVVis_ECD.docx):
  a 2-3 page report template that covers Sessions 1 and 2.
- [`Student_Report_Template_VCD.docx`](Student_Report_Template_VCD.docx): a
  2-3 page report template for Session 3.

## How the labs connect to the workflow

Every session uses the same Q-Chem pipeline scripts as the rest of the
repository: `1-qchem-init-opt.sh` through `6-qchem-tddft.sh` (electronic
spectroscopy) or `6-qchem-vcd.sh` (vibrational spectroscopy), with conformer
search (`2-qchem-conf-search.sh` or `2b-qchem-crest-conf-search.sh`),
splitting (`3-qchem-conf-split.sh` or `3b-qchem-crest-conf-split.sh`),
solvent-phase re-optimization (`4-qchem-solvent-opt.sh`), and Boltzmann
weighting (`5-qchem-boltzmann-weight.sh`). No additional dependencies are
introduced beyond what the workflow already requires.

See the [main README](../README.md) for installation, the full flag
reference, supported functionals, basis sets, and solvents, and the
molecule-list format the scripts expect.

## Editing for your course

The handouts and templates are Microsoft Word `.docx` files so that
instructors can revise wording, exercises, tables, and figures without
running a build step. Adapt the molecules, retarget for a different
audience, or modify the discussion questions as your course requires.

## License

All teaching materials are released under the MIT License along with the
rest of the repository. See [`../LICENSE`](../LICENSE).
