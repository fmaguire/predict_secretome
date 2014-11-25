#!/usr/bin/env python3
import unittest
import os
import subprocess
import predict_secretome.utils as utils

class testAncillary(unittest.TestCase):
    """
    Unittests of ancillary functions specifically 
    which and print_verbose
    """

    def test_which(self):

        output_pass = utils.which("test_files", 'test_exec')
        self.assertIs(output_pass, True)

        output_fail = utils.which("test_file", "fake_exec")
        self.assertIs(output_fail, False)


class testDependencies(unittest.TestCase):
    """
    Unittest to check dependencies work and run as expected
    """
    def setUp(self):

        execs =  ["signalp", 
                  "tmhmm",   
                  "targetp",
                  "runWolfPsortSummary"]

        path = os.path.join('dependencies', 'bin')   

        self.dependencies = {exe: path+'/'+exe \
                                for exe in execs}

        self.test_fas = os.path.join("test", 
                                     "test_files", 
                                     "test.fas")
        os.chdir(os.path.join(".."))

    def tearDown(self):
        os.chdir('test')

    def test_check_dependencies(self):
        ret = utils.check_dependencies('dependencies/bin')
        self.assertIs(True, ret)


    def test_signalp(self):

        signalp = self.dependencies['signalp']

        expected_sigpep_removed = os.path.join('test',
                                               'test_files',
                                               'expected_sigpep_removed')
        
        sigpep_removed = os.path.join('test', 
                                      'test_files', 
                                      'sigpep_removed')

        cmd = "{0} -t euk -f short -m {1} {2}".format(signalp,
                                                      sigpep_removed,
                                                      self.test_fas)
        with open(os.devnull, 'w') as null:
            retcode = subprocess.call(cmd.split(), stdout=null)

        self.assertEqual(retcode, 0)

        with open(expected_sigpep_removed, 'r') as expected_fh:
            expected_sigpep = expected_fh.readlines()

        with open(sigpep_removed, 'r') as actual_fh:
            actual_sigpep = actual_fh.readlines()

        self.assertEqual(expected_sigpep, actual_sigpep)

        os.remove(sigpep_removed)


    def test_tmhmm(self):
        
        tmhmm = self.dependencies['tmhmm']

        expected_tmhmm_output = os.path.join('test',
                                             'test_files',
                                             'expected_tmhmm_output')

        actual_tmhmm_output = os.path.join('test', 
                                           'test_files',
                                           'actual_tmhmm_output')

        cmd = "{0} {1}".format(tmhmm, self.test_fas)

        with open(actual_tmhmm_output, 'w') as out_fh:
            retcode = subprocess.call(cmd.split(), stdout=out_fh)

        self.assertEqual(retcode, 0)


        with open(expected_tmhmm_output, 'r') as expected_fh:
            expected_tmhmm = expected_fh.readlines()

        with open(actual_tmhmm_output, 'r') as actual_fh:
            actual_tmhmm = actual_fh.readlines()

        self.assertEqual(expected_tmhmm, actual_tmhmm)

        os.remove(actual_tmhmm_output)


    def test_targetp(self):
        
        targetp = self.dependencies['targetp']

        expected_targetp_output = os.path.join('test',
                                               'test_files',
                                               'expected_targetp_output')

        actual_targetp_output = os.path.join('test', 
                                             'test_files',
                                             'actual_targetp_output')

        cmd = "{0} {1}".format(targetp, self.test_fas)

        with open(actual_targetp_output, 'w') as out_fh:
            retcode = subprocess.call(cmd.split(), stdout=out_fh)

        self.assertEqual(retcode, 0)

        with open(expected_targetp_output, 'r') as expected_fh:
            expected_target = expected_fh.readlines()

        with open(actual_targetp_output, 'r') as actual_fh:
            actual_targetp = actual_fh.readlines()


        self.assertEqual(expected_target, actual_targetp)

        os.remove(actual_targetp_output)


    def test_wolfpsort(self):
        
        wolfpsort = self.dependencies['runWolfPsortSummary']

        expected_wolfp_output = os.path.join('test',
                                             'test_files',
                                             'expected_wolfp_output')

        actual_wolfp_output = os.path.join('test', 
                                           'test_files',
                                           'actual_wolfp_output')

        cmd = "{0} fungi < {1}".format(wolfpsort, 
                                       self.test_fas)

        with open(actual_wolfp_output, 'w') as out_fh:
            retcode = subprocess.call(cmd, stdout=out_fh, shell=True)

        self.assertEqual(retcode, 0)

        with open(expected_wolfp_output, 'r') as expected_fh:
            expected_wolfp = expected_fh.readlines()

        with open(actual_wolfp_output, 'r') as actual_fh:
            actual_wolfp = actual_fh.readlines()


        self.assertEqual(expected_wolfp, actual_wolfp)

        os.remove(actual_wolfp_output)

class testFormatFasta(unittest.TestCase):

  def setUp(self):
      os.chdir('..')
      self.test_fas = os.path.join("test", 
                                   "test_files", 
                                   "test.fas")

  def tearDown(self):
      os.chdir('test')

  def test_format_fasta(self):
      expected_formatted_output = os.path.join('test',
                                               'test_files',
                                               'expected_formatted_fasta')

      tmp_dir = os.path.join('test', 'test_files')

      mappings, formatted_fasta = utils.format_fasta(self.test_fas, tmp_dir)
    

      with open(expected_formatted_output, 'r') as expected_fh:
          expected_fasta = expected_fh.readlines()

      with open(formatted_fasta, 'r') as actual_fh:
            actual_fasta = actual_fh.readlines()

      self.assertEqual(expected_fasta, actual_fasta)

      self.assertEqual(len(mappings), 10)

  def test_renaming_fasta(self):
      # test renaming fasta using dict
      self.fail('finish the test')
    
if __name__=='__main__':
    unittest.main()