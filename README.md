bd-anchors: Bidirectional String Anchors
===

Bidirectional string anchors (bd-anchors, for short) is a new string sampling mechanism. Given a positive integer <b>ℓ</b>, the mechanism selects the lexicographically smallest rotation in every sliding  window of length <b>ℓ</b> of the input text. 

Bd-anchors samples are approximately uniform, locally consistent, and computable in <b>O(n)</b> time, for any input text of length <b>n</b> and any <b>ℓ</b> --- our current implementation supports an <b>O(nℓ)</b>-time construction. 

Our experiments using several datasets show that the bd-anchors sample sizes decrease proportionally to <b>ℓ</b>; and that these sizes are competitive to or smaller than the minimizers sample sizes using the analogous sampling parameters. For instance, for the Chromosome 1 of human genome, which is of length <b>n = 230,481,390</b>, and <b>ℓ = 500</b> (resp. <b>1000</b>), the set <b>A</b> of order-<b>ℓ</b> bd-anchors is of size <b>2,385,390</b> (resp. <b>1,362,218</b>). We inject bd-anchors in two problems:

<b>Indexing for On-Line Pattern Searches</b>: Our index has size <b>n</b> bytes + <b>O(|A|)</b> words and supports <b>locate</b> operations for any pattern of length at least <b>ℓ</b> in near-optimal time --- the time supported in the current implementation is not bounded. To compile the program, change to directory <b>bd-index</b> and follow the instructions given in file INSTALL.

<b>Top-K Similarity Search under Edit Distance</b>: To compile the program, change to directory <b>bd-search</b> and follow the instructions given in file INSTALL.

When publishing work that is based on the results from bd-anchors please cite:
```
G. Loukides and S. P. Pissis:
Bidirectional String Anchors: a New String Sampling Mechanism. (in press) 
ESA 2021
```

<b>License</b>: GNU GPLv3 License; Copyright (C) 2021 Grigorios Loukides and Solon P. Pissis.
