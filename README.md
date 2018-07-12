# EditTransducer: a edit transducer library for Python

This library provides an implementation of edit transducers (using the
"two-factor" construction described
[here](http://openfst.org/twiki/bin/view/FST/FstExamples)) and two
straightforward extensions:
[Levenshtein distance](https://en.wikipedia.org/wiki/Levenshtein_distance) and
[Levenshtein automata](https://en.wikipedia.org/wiki/Levenshtein_automaton).

The library has only one dependency outside of the standard library:
[Pynini](http://pynini.opengrm.org) 2.0 or better. Unfortunately, as of writing
this is not available from PyPi and thus has to be installed manually.

For usage information, see the in-module docstrings and the unit tests.

There is an accompanying tutorial, and I'll add a link once it's ready.
