#!/usr/bin/env python3

from predict_secretome.test.base import CoreTestingClass
from predict_secretome import utils

import os

class TestUtils(CoreTestingClass):
    """
    Class for testing the utility functions
    """

    def test_which(self):
        """
        Test the "which" accessory function
        """

        output_pass = utils.which(self.test_dir,
                                  'test_exec')
        self.assertIs(output_pass, True)


        output_fail = utils.which(self.test_dir, "fake_exec")
        self.assertIs(output_fail, False)


    def test_check_dependencies(self):
        """
        Make sure the dependency checking function works
        """
        ret = utils.check_dependencies(self.dep_path)
        self.assertIs(True, ret)

