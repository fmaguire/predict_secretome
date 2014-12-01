#!/usr/bin/env python3
import unittest
import os
import subprocess
import predict_secretome.utils as utils
import warnings

class coreTestingClass(unittest.TestCase):

    def compare_files(self, actual_fpath, expected_fpath):

        with open(expected_fpath, 'r') as expected_fh:
           expected_output = expected_fh.readlines()

        with open(actual_fpath, 'r') as actual_fh:
           actual_output = actual_fh.readlines()

        self.assertEqual(expected_output, actual_output)

        os.remove(actual_fpath)

    def compared_unordered_files(self, actual_fpath, expected_fpath):
        with open(expected_fpath, 'r') as expected_fh:
            expected_output = set(acc.strip() for acc in expected_fh.readlines())

        with open(actual_fpath, 'r') as actual_fh:
            actual_output = set(acc.strip() for acc in actual_fh.readlines())

        self.assertEqual(expected_output, actual_output)

        os.remove(actual_fpath)


    def compare_returned_list_to_output(self, list_of_output, output_fn):
        with open(output_fn, 'r') as out_fh:
           output = set(acc.strip() for acc in out_fh.readlines())

        self.assertEqual(set(list_of_output), output)


class testDependencies(coreTestingClass):
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


    def test_which(self):

        output_pass = utils.which("test/test_files", 'test_exec')
        self.assertIs(output_pass, True)


        output_fail = utils.which("test/test_files", "fake_exec")
        self.assertIs(output_fail, False)


    def test_check_dependencies(self):
        ret = utils.check_dependencies('dependencies/bin')
        self.assertIs(True, ret)


    def run_dependencies_and_check_output(self, expected, actual, cmd_str, dependency):

        binary = self.dependencies[dependency]

        expected_fn = os.path.join('test',
                                   'test_files',
                                    expected)

        actual_fn = os.path.join('test',
                                 'test_files',
                                  actual)

        cmd = cmd_str.format(binary, self.test_fas)

        with open(actual_fn, 'w') as out_fh:
            retcode = subprocess.call(cmd, stdout=out_fh, shell=True)


        self.assertEqual(retcode, 0)

        self.compare_files(actual_fn, expected_fn)


    def test_tmhmm(self):

        self.run_dependencies_and_check_output('expected_tmhmm_output',
                                               'actual_tmhmm_output',
                                               '{0} {1}',
                                               'tmhmm')


    def test_targetp(self):

        self.run_dependencies_and_check_output('expected_targetp_output',
                                               'actual_targetp_output',
                                               '{0} {1}',
                                               'targetp')


    def test_wolfpsort(self):

        self.run_dependencies_and_check_output('expected_wolfp_output',
                                               'actual_wolfp_output',
                                               '{0} fungi < {1}',
                                               'runWolfPsortSummary')


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

        self.compare_files(sigpep_removed, expected_sigpep_removed)



class testFormatFasta(coreTestingClass):

  def setUp(self):
      self.test_fas = os.path.join("test",
                                   "test_files",
                                   "test.fas")

      self.expected_formatted_fasta = os.path.join('test',
                                                   'test_files',
                                                    'expected_formatted_fasta')

  def test_format_fasta(self):

      tmp_dir = os.path.join('test', 'test_files')

      mappings, formatted_fasta = utils.format_fasta(self.test_fas, tmp_dir)

      self.compare_files(formatted_fasta, self.expected_formatted_fasta)

      self.assertEqual(len(mappings), 10)


  def test_renaming_fasta(self):
      # test renaming fasta using dict

      with open(self.expected_formatted_fasta) as formatted_fh:
          accessions = [acc.rstrip('\n').lstrip('>') for \
                           acc in formatted_fh.readlines() if acc.startswith('>')]

      mappings = {'1c3c35783faef7c620b': 'R1qual32.paired_(paired)_contig_1567',
                  '8281d5113ad095094d4': 'PM50_trimmed_paired_qual32_contig_18412',
                  'df82f90060f7e0d6dba': 'PM30_NoIndex_L003_R1_001_(paired)_trimmed_(paired)_contig_2634',
                  'faa0e74320c2db42758': 'PM30_NoIndex_L003_R1_001_(paired)_trimmed_(paired)_contig_8670',
                  'b4afdb58cd630e0b6a9': 'PM50_trimmed_paired_qual32_contig_12144',
                  '747588efe930ab229a5': 'PM50_trimmed_paired_qual32_contig_10408',
                  '726c9585f9462b87b5f': 'PM50_trimmed_paired_qual32_contig_5857',
                  'ae875b7838abc202013': 'PM30_NoIndex_L003_R1_001_(paired)_trimmed_(paired)_contig_72552',
                  '83440bd1de03a46193e': 'PM50_trimmed_paired_qual32_contig_16575',
                  '21cbfe90785a9f098f8': 'R1qual32.paired_(paired)_contig_2534'}

      utils.generate_output(self.expected_formatted_fasta, accessions, mappings, 'test')

      output_fn = 'test_conservative_predicted_secretome.fasta'

      self.compare_files(output_fn, self.test_fas)


class testSecretome(coreTestingClass):

    def setUp(self):
        self.list_1 = ['A', 'B']
        self.list_2 = ['A', 'B', 'C']
        self.list_3 = ['A', 'B', 'D', 'E']
        self.list_4 = ['A', 'B', 'E', 'F']

        self.tmp_dir = 'test_files'


    def test_conservative(self):
        expected_permissive_list_fn = os.path.join(self.tmp_dir,
                                                   "expected_conservative_secretome_accessions.txt")

        accession_list = utils.secretome(self.list_1,
                                         self.list_2,
                                         self.list_3,
                                         self.list_4,
                                         self.tmp_dir)

        self.assertEqual({'A', 'B'}, accession_list)

    def test_permissive(self):
        expected_permissive_list_fn = os.path.join(self.tmp_dir,
                                                   "expected_permissive_secretome_accessions.txt")

        accession_list = utils.secretome(self.list_1,
                                         self.list_2,
                                         self.list_3,
                                         self.list_4,
                                         self.tmp_dir,
                                         conservative=False)

        self.assertEqual({'A', 'B', 'C', 'D', 'E', 'F'}, accession_list)


    def test_null_output(self):
        null = []


        with warnings.catch_warnings(record=True) as w:


            accession_list = utils.secretome(self.list_1,
                                             self.list_2,
                                             self.list_3,
                                             null,
                                             self.tmp_dir)

            self.assertEqual(len(w), 1)
            self.assertIs(w[-1].category, UserWarning)
            self.assertEqual(str(w[-1].message),
                            "No secreted proteins found using {0} setting".format('conservative'))

            self.assertEqual(set(), accession_list)


class testDependencyParsing(coreTestingClass):

    def setUp(self):
        self.formatted_fasta = os.path.join('test',
                                            'test_files',
                                            'expected_formatted_fasta')
        self.tmp_dir = os.path.join('test', 'test_files')

        self.test_fas = os.path.join(self.tmp_dir, 'test.fas')


        #self.expected_seqs_with_sigpep_removed = os.path.join(self.tmp_dir, '')



    def test_tmhmm_func(self):
        expected_acc_without_tm_domains = os.path.join(self.tmp_dir,
                                                      'test_without_tm_domains.txt')

        actual_acc_without_tm_domains = utils.tmhmm(self.test_fas, '')

        self.assertEqual(["PM50_trimmed_paired_qual32_contig_10408",
                          "R1qual32.paired_(paired)_contig_1567",
                          "PM50_trimmed_paired_qual32_contig_16575",
                          "PM50_trimmed_paired_qual32_contig_18412",
                          "PM30_NoIndex_L003_R1_001_(paired)_trimmed_(paired)_contig_72552",
                          "PM30_NoIndex_L003_R1_001_(paired)_trimmed_(paired)_contig_2634",
                          "R1qual32.paired_(paired)_contig_2534",
                          "PM30_NoIndex_L003_R1_001_(paired)_trimmed_(paired)_contig_8670"],
                          actual_acc_without_tm_domains)


    def test_signalp_func(self):
        expected_accessions_with_sig_pep = ['PM50_trimmed_paired_qual32_contig_10408',
                                            'R1qual32.paired_(paired)_contig_1567',
                                            'PM50_trimmed_paired_qual32_contig_16575',
                                            'PM50_trimmed_paired_qual32_contig_18412',
                                            'PM50_trimmed_paired_qual32_contig_12144',
                                            'R1qual32.paired_(paired)_contig_2534']

        expected_mature_seq_fp = os.path.join(self.tmp_dir, 'expected_signalp_mature_seqs.fasta')

        expected_full_sequences_with_sigpep_fp = os.path.join(self.tmp_dir,
                                                              'expected_signalp_full_seqs_with_sigpep.fasta')

        actual_mature_seqs_fp, \
        actual_accessions_with_sig_pep, \
        actual_full_sequences_with_sigpep_fp = utils.signalp(self.test_fas, '')

        self.assertEqual(actual_accessions_with_sig_pep, expected_accessions_with_sig_pep)

        self.compare_files(actual_full_sequences_with_sigpep_fp, expected_full_sequences_with_sigpep_fp)

        self.compare_files(actual_mature_seqs_fp, expected_mature_seq_fp)



    def test_wolfpsort_func(self):

        expected_extracellular_accessions = os.path.join(self.tmp_dir,
                                                         'test_extr_acc.txt')

        actual_extracellular_accessions = utils.wolfpsort(self.test_fas, '')

        self.assertEqual(["PM50_trimmed_paired_qual32_contig_10408",
                         "R1qual32.paired_(paired)_contig_1567"],
                          actual_extracellular_accessions)


    def test_targetp_func(self):

        expected_secreted_accessions = os.path.join(self.tmp_dir,
                                                         'test_s_acc_truncated.txt')

        actual_secreted_accessions = utils.targetp(self.test_fas, '')

        self.assertEqual(["PM50_trimmed_paired_",
                          "R1qual32.paired_(pai",
                          "PM50_trimmed_paired_",
                          "PM50_trimmed_paired_",
                          "PM30_NoIndex_L003_R1",
                          "R1qual32.paired_(pai",
                          "PM30_NoIndex_L003_R1"],
                          actual_secreted_accessions)


if __name__=='__main__':
    unittest.main()
