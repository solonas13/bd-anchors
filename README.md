bd-anchors: Bidirectional String Anchors
===

Bidirectional string anchors (bd-anchors, for short) is a new string sampling mechanism. Given a positive integer <b>ℓ</b>, our mechanism selects the lexicographically smallest rotation in every sliding  window of length <b>ℓ</b> of the input text. 

Bd-anchors samples are approximately uniform, locally consistent, and computable in linear time --- our current implementation supports <b>O(nℓ)</b>-time computation, where <b>n</b> is the length of the input text. 

Our experiments using several datasets show that the bd-anchors sample sizes decrease proportionally to <b>ℓ</b>; and that these sizes are competitive to or smaller than the minimizers sample sizes using the analogous sampling parameters. For instance, for the Chromosome 1 of human genome, which is of length <b>n = 230,481,390</b>, and <b>ℓ = 500</b>, the set <b>A</b> of bd-anchors is of size <b>2,385,390</b>.

<b>Index</b>: the current implementation constructs an index of size <b>n</b> bytes + <b>O(|A|)</b> 64-bit words. It supports <b>locate</b> operations for any pattern of length at least <b>ℓ</b>. To compile the program, please change to directory <b>index</b> and follow the instructions given in file INSTALL.

When publishing work that is based on the results from bd-anchors please cite:
```
G. Loukides and S. P. Pissis:
Bidirectional String Anchors: a New String Sampling Mechanism. 
ESA 2021
```

<b>License</b>: GNU GPLv3 License; Copyright (C) 2021 Grigorios Loukides and Solon P. Pissis.
