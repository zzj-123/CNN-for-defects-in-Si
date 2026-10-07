# CNN for atomic-defect classification in silicon

Historical research by **Zijian Zhou (2019-2021)**. This repository is being documented in October 2026 using existing reports, slides and code.

## Research summary

**[Read the illustrated research summary (PDF)](docs/research-summary.pdf)**

The project investigated neural-network classification of atomic-defect conditions from simulated electronic band-structure data. Alongside network training and evaluation, I explored what input information supported classification and how model behavior changed across conditions.

The summary presents three connected analyses:

- **Input occlusion:** regional masking to probe CNN responses, including differences in sensitivity patterns between models.
- **Representation comparisons:** histogram, shuffled and thresholded inputs to investigate discriminative information and distinguish task-relevant information from the mechanisms used by a particular model.
- **Generalization:** comparisons across simulation conditions, spectral resolution and energy windows.

My contributions included dataset generation and preprocessing, neural-network implementation and training, evaluation, and the exploratory analyses described above. The project used SIESTA, LAMMPS and Python, with fully connected networks and CNNs, including PyTorch CNN implementations.

## Historical scope

The results and figures come from the original project. The October 2026 update organizes this material; it does not introduce new experiments. Input-level analyses did not establish internal circuits or resolve the original CNN's mechanisms. Successful classification using reduced input information was not treated as proof that the original CNN used the same information.

## Repository contents

- `docs/research-summary.pdf`: the illustrated historical research summary.
- `figures/`: original figure assets used in the summary.
- `docs/figure-provenance.json`: source slide and image locators for the reproduced figures.
- `rdf.py`: the original 2020 radial-distribution preprocessing script, retained unchanged. It expects project-specific LAMMPS inputs.
- `test.txt`: an existing empty historical file, retained unchanged.

The full simulation datasets, model checkpoints and training environment are not distributed in this documentation update.
