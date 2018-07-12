# Copyright (C) 2017 Kyle Gorman
#
# Permission is hereby granted, free of charge, to any person obtaining a
# copy of this software and associated documentation files (the "Software"),
# to deal in the Software without restriction, including without limitation the
# rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
# sell copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

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
"""

from __future__ import division

from pynini import compose
from pynini import invert
from pynini import NO_STATE_ID
from pynini import shortestdistance
from pynini import shortestpath
from pynini import string_map
from pynini import transducer
from pynini import union


DEFAULT_INSERT_COST = 1
DEFAULT_DELETE_COST = 1
DEFAULT_SUBSTITUTE_COST = 1


class LatticeError(Exception):
  pass


class EditTransducer(object):
  """Factored edit transducer.

  This class stores the two factors of an finite-alphabet edit transducer and
  supports insertion, deletion, and substitution operations with user-specified
  costs.

  Note that the cost of substitution must be less than the cost of insertion
  plus the cost of deletion or no optimal path will include substitution.
  """

  # Reserved labels for edit operations.
  DELETE = "<delete>"
  INSERT = "<insert>"
  SUBSTITUTE = "<substitute>"

  def __init__(self,
               alphabet,
               insert_cost=DEFAULT_INSERT_COST,
               delete_cost=DEFAULT_DELETE_COST,
               substitute_cost=DEFAULT_SUBSTITUTE_COST):
    """Constructor.

    Args:
      alphabet: edit alphabet (an iterable of strings).
      insert_cost: the cost for the insertion operation.
      delete_cost: the cost for the deletion operation.
      substitute_cost: the cost for the substitution operation.
    """
    # Left factor; note that we divide the edit costs by two because they also
    # will be incurred when traversing the right factor.
    match = union(*alphabet).optimize(True)
    i_insert = transducer("", "[{}]".format(self.INSERT),
                          weight=insert_cost / 2).optimize(True)
    i_delete = transducer(match, "[{}]".format(self.DELETE),
                          weight=delete_cost / 2).optimize(True)
    i_substitute = transducer(match, "[{}]".format(self.SUBSTITUTE),
                              weight=substitute_cost / 2).optimize(True)
    i_ops = union(match, i_insert, i_delete, i_substitute).optimize(True)
    # Right factor; this is constructed by inverting the left factor (i.e.,
    # swapping the input and output labels), then swapping the insert and delete
    # labels on what is now the input side.
    o_ops = invert(i_ops)
    syms = o_ops.input_symbols()
    insert_label = syms.find(self.INSERT)
    delete_label = syms.find(self.DELETE)
    o_ops.relabel_pairs(ipairs=((insert_label, delete_label),
                                (delete_label, insert_label)))
    # Computes the closure for both sets of ops.
    self._e_i = i_ops.closure().optimize(True)
    self._e_o = o_ops.closure().optimize(True)
 
  @staticmethod
  def check_wellformed_lattice(lattice):
    """Raises an error if the lattice is empty.

    Args:
      lattice: A lattice FST.

    Raises:
      LatticeError: Lattice is empty.
    """
    if lattice.start() == NO_STATE_ID:
      raise LatticeError("Lattice is empty")

  def _create_lattice(self, iset, oset):
    """Creates edit lattice for a pair of input/output strings or acceptors.

    Args:
      iset: input string or acceptor
      oset: output string or acceptor.

    Returns:
      A lattice FST.
    """
    l_i = compose(iset, self._e_i)
    l_o = compose(self._e_o, oset)
    lattice = compose(l_i, l_o)
    EditTransducer.check_wellformed_lattice(lattice)
    return lattice


class LevenshteinDistance(EditTransducer):
  """Edit transducer augmented with a distance calculator."""

  def distance(self, iset, oset):
    """Computes minimum distance.

    This method computes, for a pair of input/output strings or acceptors, the
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
               delete_cost=DEFAULT_DELETE_COST,
               substitute_cost=DEFAULT_SUBSTITUTE_COST):
    super(LevenshteinAutomaton, self).__init__(alphabet, insert_cost,
                                               delete_cost, substitute_cost)
    # Compiles lexicon and composes the right factor with it.
    compiled_lexicon = string_map(lexicon)
    self._l_o = compose(self._e_o, compiled_lexicon)
    self._l_o.optimize(True)

  def _create_levenshtein_automaton_lattice(self, query):
    """Constructs a lattice for a query string.

    Args:
      query: input string or acceptor.

    Returns:
      A lattice FST.
    """
    l_i = compose(query, self._e_i)
    lattice = compose(l_i, self._l_o)
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
    lattice = self._create_levenshtein_automaton_lattice(query)
    # For implementation reasons, the shortest path (when k = 1) is in reverse
    # state order, so we perform a topological sort ahead of time.
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
      An iterator of the closest strings in the lexicon.
    """
    lattice = self._create_levenshtein_automaton_lattice(query)
    # Prunes all paths whose weights are worse than the best path.
    lattice.prune(weight=0).project(True).optimize(True)
    return lattice.paths().ostrings()
