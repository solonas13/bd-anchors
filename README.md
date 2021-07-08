bd-anchors: Bidirectional String Anchors
===

Bidirectional string anchors (bd-anchors, for short) is a new string sampling mechanism. Given a positive integer <b>ell</b>, our mechanism selects the lexicographically smallest rotation in every sliding  window of length <b>ell</b>. 

Bd-anchors samples are approximately uniform, locally consistent, and computable in linear time --- our current implementation supports <b>O(n.ell)</b>-time computation, where <b>n</b> is the length of the input text. 

Our experiments using several datasets demonstrate that the bd-anchors sample sizes decrease proportionally to <b>ell</b>; and that these sizes are competitive to or smaller than the minimizers sample sizes using the analogous sampling parameters.

<b>Index</b>: the current implementations constructs an index of size <b>n</b> bytes + <b>O(|A|)</b> words, where <b>A</b> is the set of bd-anchors. To compile the program, please change to directory <b>index</b> and follow the instructions given in file INSTALL.

When publishing work that is based on the results from bd-anchors please cite:
```
G. Loukides and S. P. Pissis:
Bidirectional String Anchors: a New String Sampling Mechanism. 
ESA 2021
```

<b>License</b>: GNU GPLv3 License; Copyright (C) 2014 Grigorios Loukides and Solon P. Pissis.
