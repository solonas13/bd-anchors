bd-anchors: Bidirectional String Anchors
===

Bidirectional string anchors (bd-anchors) is a new string sampling mechanism. Given a positive integer <b>ℓ</b>, the mechanism selects the leftmost lexicographically smallest rotation in every sliding  window of length <b>ℓ</b> of the input text. 

Bd-anchors samples are approximately uniform, locally consistent, and computable in <b>O(n)</b> time, for any input text of length <b>n</b> and any <b>ℓ</b> --- our current implementation supports an <b>O(nℓ)</b>-time construction. 

Our experiments using several datasets show that the bd-anchors sample sizes decrease proportionally to <b>ℓ</b>; and that these sizes are competitive to or smaller than the minimizers sample sizes using the analogous sampling parameters. For instance, for the Chromosome 1 of human genome, which is of length <b>n = 230,481,390</b>, and <b>ℓ = 500</b> (resp. <b>1000</b>), the set <b>A</b> of order-<b>ℓ</b> bd-anchors is of size <b>1,560,882</b> (resp. <b>897,953</b>). 

<b>Constructing the Sample</b>: Our current implementation takes <b>O(nℓ)</b> time. To compile the program, change to directory <b>bd-construct</b> and follow the instructions given in file INSTALL. 

We inject <b>bd-anchors</b> in two problems:

<b>Text Indexing</b>: Our index has size <b>n</b> bytes + <b>O(|A|)</b> integers and supports <b>locate</b> operations for any pattern of length at least <b>ℓ</b> in near-optimal time (<b>bd-index-grid</b>) --- the time supported in the <b>bd-index</b> implementation is not bounded but this implementation is faster when the text is not very long (e.g.~200MB). To compile the program, change to directory <b>bd-index</b> or <b>bd-index-grid</b> and follow the instructions given in file INSTALL.

<b>Top-K Similarity Search under Edit Distance</b>: To compile the program, change to directory <b>bd-search</b>  and follow the instructions given in file INSTALL.

When publishing work that is based on the results from bd-anchors please cite:
```
G. Loukides, S. P. Pissis, M. Sweering:
Bidirectional String Anchors for Improved Text Indexing and Top-K Similarity Search. 
IEEE Trans. Knowl. Data Eng. DOI: 10.1109/TKDE.2022.3231780
```
```
G. Loukides and S. P. Pissis:
Bidirectional String Anchors: a New String Sampling Mechanism. 
ESA 2021: 64:1-64:21. DOI: 10.4230/LIPIcs.ESA.2021.64
```

<b>License</b>: GNU GPLv3 License; Copyright (C) 2021 Grigorios Loukides and Solon P. Pissis.
