# Encoding: UTF-8
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

"""Tests edit transducer classes (specifically, the Levenshtein automaton)."""

import string
import unittest

import edit_transducer


class LevenshteinAutomatonTest(unittest.TestCase):

  @classmethod
  def setUpClass(cls):
    cheese_lexicon = ("tilsit", "caerphilly", "stilton", "gruyere", "emmental",
                      "liptauer", "lancashire", "cheshire", "brie", "roquefort",
                      "savoyard", "boursin", "camembert", "gouda", "edam",
                      "caithness", "wensleydale", "gorgonzola", "parmesan",
                      "mozzarella", "fynbo", "cheddar", "ilchester",
                      "limburger")
    cls.automaton = edit_transducer.LevenshteinAutomaton(
        string.ascii_lowercase, cheese_lexicon)
    cls.distance = edit_transducer.LevenshteinDistance(string.ascii_lowercase)

  def query_and_distance(self, query, expected_closest, expected_distance):
    closest = self.automaton.closest_match(query)
    self.assertEqual(expected_closest, closest)
    distance = self.distance.distance(query, closest)
    self.assertEqual(expected_distance, distance)

  ## Tests using query_and_distance helper.

  def testMatch(self):
    self.query_and_distance("stilton", "stilton", 0.0)

  def testInsertion(self):
    self.query_and_distance("mozarela", "mozzarella", 2.0)

  def testDeletion(self):
    self.query_and_distance("emmenthal", "emmental", 1.0)

  def testSubstitution(self):
    self.query_and_distance("bourzin", "boursin", 1.0)

  def testMixedEdit(self):
    self.query_and_distance("rockford", "roquefort", 4.0)

  ## Other tests

  def testClosestMatchFindsExactMatch(self):
    res = self.automaton.closest_matches("cheddar")
    self.assertListEqual(["cheddar"], list(res))

  def testClosestMatchReturnsMultiple(self):
    res = self.automaton.closest_matches("cheese")
    self.assertListEqual(["cheddar", "cheshire"], list(res))

  def testOutOfAlphabetQueryRaisesError(self):
    with self.assertRaises(edit_transducer.LatticeError):
      unused_closest = self.automaton.closest_match("Gruyère")


if __name__ == "__main__":
  unittest.main()
