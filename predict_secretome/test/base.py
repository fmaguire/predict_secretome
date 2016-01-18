#!/usr/bin/env python3

import unittest
import os

class CoreTestingClass(unittest.TestCase):
    """
    Main superclass with utility functions
    """
    dep_path = os.path.join('predict_secretome', 'dependencies', 'bin')
    test_dir = os.path.join('predict_secretome', 'test', 'test_files')


    def compare_files(self, actual_fpath, expected_fpath):
        """
        Test that two files are the same
        """

        with open(expected_fpath, 'r') as expected_fh:
            expected_output = expected_fh.readlines()

        with open(actual_fpath, 'r') as actual_fh:
            actual_output = actual_fh.readlines()

        self.assertEqual(expected_output, actual_output)

        os.remove(actual_fpath)


    def compared_unordered_files(self, actual_fpath, expected_fpath):
        """
        Test that two files are the same when order in the file doesn't
        matter
        """
        with open(expected_fpath, 'r') as expected_fh:
            expected_output = set(acc.strip() for acc in expected_fh.readlines())

        with open(actual_fpath, 'r') as actual_fh:
            actual_output = set(acc.strip() for acc in actual_fh.readlines())

        self.assertEqual(expected_output, actual_output)

        os.remove(actual_fpath)


    def compare_returned_list_to_output(self, list_of_output, output_fn):
        """
        Test that a list and a file have the same contents
        """
        with open(output_fn, 'r') as out_fh:
            output = set(acc.strip() for acc in out_fh.readlines())

        self.assertEqual(set(list_of_output), output)

