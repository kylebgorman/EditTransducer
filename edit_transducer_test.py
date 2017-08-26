# Encoding: UTF-8
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

  def query_and_distance(self, query, expected_closest, expected_distance):
    closest = self.automaton.closest_match(query)
    self.assertEqual(expected_closest, closest)
    distance = self.automaton.distance(query, closest)
    self.assertEqual(expected_distance, distance)

  def testMatch(self):
    self.query_and_distance("stilton", "stilton", 0.0)

  def testInsertion(self):
    self.query_and_distance("cheeshire", "cheshire", 1.0)

  def testDeletion(self):
    self.query_and_distance("mozarela", "mozzarella", 2.0)

  def testMixedEdit(self):
    self.query_and_distance("rockford", "roquefort", 7.0)

  def testOutOfAlphabetQueryRaisesError(self):
    with self.assertRaises(edit_transducer.Error):
      unused_closest = self.automaton.closest_match("Gruyère")


if __name__ == "__main__":
  unittest.main()
