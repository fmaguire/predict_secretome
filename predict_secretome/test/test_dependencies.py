#!/usr/bin/env python3

from predict_secretome.test.base import CoreTestingClass

from predict_secretome import utils

import os
import subprocess
import shutil

class TestDependencies(CoreTestingClass):
    """
    Unittest to check dependencies work and run as expected
    """
    def setUp(self):

        execs = ["signalp",
                 "tmhmm",
                 "targetp",
                 "runWolfPsortSummary"]


        self.dependencies = {exe: self.dep_path+'/'+exe \
                                for exe in execs}


        self.test_fas = os.path.join(self.test_dir,
                                    "test.fas")


    def run_dependencies_and_check_output(self, expected, actual, cmd_str, dependency):
        """
        Utility function which runs the dependency supplied and checks
        the output
        """

        binary = self.dependencies[dependency]

        expected_fn = os.path.join(self.test_dir,
                                   expected)

        actual_fn = os.path.join(self.test_dir,
                                 actual)

        cmd = cmd_str.format(binary, self.test_fas)

        with open(actual_fn, 'w') as out_fh:
            retcode = subprocess.call(cmd, stdout=out_fh, shell=True)


        self.assertEqual(retcode, 0)

        self.compare_files(actual_fn, expected_fn)



    def test_tmhmm(self):
        """
        Test TMHMM depedency is installed and works as expected
        """

        self.run_dependencies_and_check_output('expected_tmhmm_output',
                                               'actual_tmhmm_output',
                                               '{0} {1}',
                                               'tmhmm')


    def test_targetp(self):
        """
        Test the targetp dependency is installed and works as expected
        """

        username = os.getlogin()
        home_targetp = '/home/{0}/targetp-1.1'.format(username)
        shutil.copytree('predict_secretome/dependencies/targetp-1.1', home_targetp)

        self.run_dependencies_and_check_output('expected_targetp_output',
                                               'actual_targetp_output',
                                               '{0} {1}',
                                               'targetp')
        shutil.rmtree(home_targetp)


    def test_wolfpsort(self):
        """
        Test that the wolfpsort dependency is installed and works
        """

        self.run_dependencies_and_check_output('expected_wolfp_output',
                                               'actual_wolfp_output',
                                               '{0} fungi < {1}',
                                               'runWolfPsortSummary')


    def test_signalp(self):
        """
        Test that the signalp dependency is installed and works
        """

        signalp = self.dependencies['signalp']

        expected_sigpep_removed = os.path.join(self.test_dir,
                                               'expected_sigpep_removed')

        sigpep_removed = os.path.join(self.test_dir,
                                      'sigpep_removed')

        cmd = "{0} -t euk -f short -m {1} {2}".format(signalp,
                                                      sigpep_removed,
                                                      self.test_fas)
        with open(os.devnull, 'w') as null:
            retcode = subprocess.call(cmd.split(), stdout=null)

        self.assertEqual(retcode, 0)

        self.compare_files(sigpep_removed, expected_sigpep_removed)


