# Encoding: UTF-8
"""Edit transducer classes.

Edit transducers are abstract machines used to efficiently compute edit
distance and approximate string matches.

Here, we provide three concrete classes:

* EditTransducer(object): Constructs the transducer from an input alphabet and
  cost matrix. Provides a protected `_make_lattice` method for lattice
  construction, which may be overridden by derived classes.
* LevenshteinDistance(EditTransducer): Also adds a method for computing edit distance
  from the lattice.
* LevenshteinAutomaton(LevenshteinDistance): Uses the edit transducer and an input
  vocabulary to construct a right-factored lexicon, from which one can compute
  the closest matches.

One current limitation is that the library only supports "indel" (insertion
and deletion) operations, so apparent substitutions are modeled as an insertion
followed by a deletion (or vis versa). It should be straightforward to extend
this model to support substitutions, however.
"""

from __future__ import division

from pynini import acceptor, string_map, transducer, NO_STATE_ID
from pynini import compose, invert, shortestdistance, shortestpath, union


DEFAULT_INSERT_COST = 1.
DEFAULT_DELETE_COST = 1.


class Error(Exception):
  """Generic error for this module."""
  pass


class EditTransducer(object):
  """Factored edit transducer supporting indel operations.

  The two-factor construction used here is based on the "Edit Distance"
  OpenFst exercise (http://www.openfst.org/twiki/bin/view/FST/FstExamples),
  which is much faster than the naïve (i.e., unfactored) construction.
  """

  # Reserved labels for edit operations.
  DELETE = "<delete>"
  INSERT = "<insert>"

  def __init__(self,
               alphabet,
               insert_cost=DEFAULT_INSERT_COST,
               delete_cost=DEFAULT_DELETE_COST):
    """Constructor.

    Args:
      alphabet: edit alphabet (an iterable of strings).
      insert_cost: the cost for the insertion operation.
      delete_cost: the cost for the deletion operation.
    """
    # Left factor.
    match = union(*alphabet).optimize(True)
    insert = transducer("", "[{}]".format(self.INSERT), weight=insert_cost / 2.)
    delete = transducer(match, "[{}]".format(self.DELETE),
                        weight=delete_cost / 2.).optimize(True)
    self._left_factor = union(match, insert, delete).optimize(True)
    self._left_factor.closure().optimize(True)
    # Right factor.
    self._right_factor = invert(self._left_factor)
    syms = self._right_factor.input_symbols()
    insert_label = syms.find(self.INSERT)
    delete_label = syms.find(self.DELETE)
    self._right_factor.relabel_pairs(ipairs=((insert_label, delete_label),
                                             (delete_label, insert_label)))
    self._right_factor.closure().optimize(True)

  @staticmethod
  def check_wellformed_lattice(lattice):
    """Raises an error if the lattice is empty.

    Args:
      lattice: A lattice FST.

    Raises:
      Error: Lattice is empty.
    """
    if lattice.start() == NO_STATE_ID:
      raise Error("Lattice is empty")

  def _create_lattice(self, iset, oset):
    """Creates edit lattice for a pair of input/output strings or acceptors.

    Args:
      iset: input string or acceptor
      oset: output string or acceptor.

    Returns:
      A lattice FST.
    """
    left = compose(iset, self._left_factor)
    right = compose(self._right_factor, oset)
    lattice = compose(left, right)
    EditTransducer.check_wellformed_lattice(lattice)
    return lattice


class LevenshteinDistance(EditTransducer):
  """Edit transducer augmented with a distance calculator."""

  def distance(self, iset, oset):
    """Computes minimum distance.

    This method comptues, for a pair of input/output strings or acceptors, the
    minimum edit distance according to the underlying edit transducer.

    Args:
      iset: input string or acceptor.
      oset: output string or acceptor.

    Returns:
      Minimum edit distance according to the edit transducer.
    """
    lattice = self._create_lattice(iset, oset)
    # The shortest cost from all final states to the start state is
    # equivalent to the cost of the shortest path.
    return float(shortestdistance(lattice, reverse=True)[0])


class LevenshteinAutomaton(LevenshteinDistance):
  """Edit transducer with a fixed output lexicon and closest-match methods."""

  def __init__(self,
               alphabet,
               lexicon,
               insert_cost=DEFAULT_INSERT_COST,
               delete_cost=DEFAULT_DELETE_COST):
    super(LevenshteinAutomaton, self).__init__(alphabet, insert_cost,
                                               delete_cost)
    compiled_lexicon = string_map(lexicon).optimize(True)
    self._factored_lexicon = compose(self._right_factor, compiled_lexicon)
    self._factored_lexicon.optimize(True)

  def _make_lattice(self, query):
    """Constructs a lattice for a query string.

    Args:
      query: input string or acceptor.

    Returns:
      A lattice FST.
    """
    left = compose(query, self._left_factor)
    lattice = compose(left, self._factored_lexicon)
    EditTransducer.check_wellformed_lattice(lattice)
    return lattice

  def closest_match(self, query):
    """Returns the closest string to the query in the lexicon.

    This method computes, for an input string or acceptor, the closest string
    in the lexicon according to the underlying edit transducer. In the case of
    a tie (i.e., where there are multiple closest strings), only one will be
    returned; tie breaking is deterministic but difficult to reason about and
    thus should be considered unspecified.) The `closest_matches` method can be
    used to enumerate all the ties.

    Args:
      query: input string or acceptor.

    Returns:
      The closest string in the lexicon.
    """
    lattice = self._make_lattice(query)
    # For implementational reasons, the shortest path (when k = 1) is in
    # reverse state order, so we perform a topological sort ahead of time.
    return shortestpath(lattice).topsort().stringify()

  def closest_matches(self, query):
    """Returns all of the closest strings to the query in the lexicon.

    This method returns, for an input string or acceptor, the closest strings
    in the lexicon according to the underlying edit transducer. A string is
    "closest" if it has the same edit distance as the closest string. The order
    in which the strings are returned is deterministic but difficult to reason
    about and thus should be considered unspecified.

    Args:
      query: input string or acceptor.

    Returns:
      A tuple of the closest strings in the lexicon.
    """
    lattice = self._make_lattice(query)
    # Prunes all paths whose weights are worse than the best path.
    lattice.prune(weight=0.).project(True).optimize(True)
    return tuple(lattice.paths().iter_ostrings())
