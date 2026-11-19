# External 7 Public Data Reconstruction

This notebook reconstructs the public versions of the external benchmark
datasets used in the HSA evaluation.

Seven publications are covered:

1. Schukken
2. Abt-Wang
3. Wilson
4. Nguyen
5. Shan
6. Herranz-Ors
7. Samelson

The workflow downloads public source data, normalizes gene identifiers,
builds AssayLoop-compatible screens, and verifies reconstructed gene/hit
counts against the expected benchmark definitions.

This is a public-data reconstruction and is not necessarily identical to
Odysseus's frozen external assets where those assets are unavailable.

Known caveat:
Schukken's public source permits more than one reasonable target
reconstruction. Exact reproduction of the frozen E224 194-hit target
requires the original frozen target identities.
