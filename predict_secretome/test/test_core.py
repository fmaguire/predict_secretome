#!/usr/bin/env python3
"""
Unittests for secretome prediction script/module
"""
from __future__ import print_function
import unittest
import os
import shutil
import subprocess
import predict_secretome.utils as utils
import warnings
from predict_secretome.test.base import CoreTestingClass

#class testPredictSecretomeClass(CoreTestingClass):
#    """
#    Tests for the main predictSecretome runner class
#    """
#
#    def setUp(self):
#        pass
#
#    def test_init(self):
#        pass
#




#class TestSecretome(CoreTestingClass):
#    """
#    Test the functions that combine predictions
#    """
#
#    def setUp(self):
#        """
#        Generate lists shared by tests
#        """
#        self.list_1 = ['A', 'B']
#        self.list_2 = ['A', 'B', 'C']
#        self.list_3 = ['A', 'B', 'D', 'E']
#        self.list_4 = ['A', 'B', 'E', 'F']
#
#        self.tmp_dir = 'test_files'
#
#
#    def test_conservative(self):
#        """
#        Test conservative prediction generation
#        """
#
#        accession_list = utils.secretome(self.list_1,
#                                         self.list_2,
#                                         self.list_3,
#                                         self.list_4)
#
#        self.assertEqual({'A', 'B'}, accession_list)
#
#    def test_permissive(self):
#        """
#        Test permissive prediction generation
#        """
#        accession_list = utils.secretome(self.list_1,
#                                         self.list_2,
#                                         self.list_3,
#                                         self.list_4,
#                                         conservative=False)
#
#        self.assertEqual({'A', 'B', 'C', 'D', 'E', 'F'}, accession_list)
#
#
#    def test_null_output(self):
#        """
#        Test nothing is output correctly when conservative criteria aren't
#        matched
#        """
#
#        null = []
#
#        with warnings.catch_warnings(record=True) as warn:
#            accession_list = utils.secretome(self.list_1,
#                                             self.list_2,
#                                             self.list_3,
#                                             null)
#
#            self.assertEqual(len(warn), 1)
#            self.assertIs(warn[-1].category, UserWarning)
#            self.assertEqual(str(warn[-1].message),
#                             "No secreted proteins found using {0} setting".format('conservative'))
#
#            self.assertEqual(set(), accession_list)
#
#
#class TestDependencyParsing(CoreTestingClass):
#    """
#    Test the function that call and parse the various dependencies
#    """
#
#    def setUp(self):
#        """
#        Generate shared inputs
#        """
#        self.formatted_fasta = os.path.join('test',
#                                            'test_files',
#                                            'expected_formatted_fasta')
#        self.tmp_dir = os.path.join('test', 'test_files')
#
#        self.test_fas = os.path.join(self.tmp_dir, 'test.fas')
#
#
#    def test_tmhmm_func(self):
#        """
#        Test the tmhmm function works and parses output correctly
#        """
#
#        actual_acc_without_tm_domains = utils.tmhmm(self.test_fas)
#
#        self.assertEqual(["PM50_trimmed_paired_qual32_contig_10408",
#                          "R1qual32.paired_(paired)_contig_1567",
#                          "PM50_trimmed_paired_qual32_contig_16575",
#                          "PM50_trimmed_paired_qual32_contig_18412",
#                          "PM30_NoIndex_L003_R1_001_(paired)_trimmed_(paired)_contig_72552",
#                          "PM30_NoIndex_L003_R1_001_(paired)_trimmed_(paired)_contig_2634",
#                          "R1qual32.paired_(paired)_contig_2534",
#                          "PM30_NoIndex_L003_R1_001_(paired)_trimmed_(paired)_contig_8670"],
#                         actual_acc_without_tm_domains)
#
#
#    def test_signalp_func(self):
#        """
#        Test that the signalp function correctly parses signalp output
#        """
#
#        expected_accessions_with_sig_pep = ['PM50_trimmed_paired_qual32_contig_10408',
#                                            'R1qual32.paired_(paired)_contig_1567',
#                                            'PM50_trimmed_paired_qual32_contig_16575',
#                                            'PM50_trimmed_paired_qual32_contig_18412',
#                                            'PM50_trimmed_paired_qual32_contig_12144',
#                                            'R1qual32.paired_(paired)_contig_2534']
#
#        expected_mature_seq_fp = os.path.join(self.tmp_dir, 'expected_signalp_mature_seqs.fasta')
#
#        expected_full_sequences_with_sigpep_fp = os.path.join(self.tmp_dir,
#                                                              'expected_signalp_full_seqs_with_sigpep.fasta')
#
#        actual_mature_seqs_fp, \
#        actual_accessions_with_sig_pep, \
#        actual_full_sequences_with_sigpep_fp = utils.signalp(self.test_fas, '', '', '')
#
#        self.assertEqual(actual_accessions_with_sig_pep, expected_accessions_with_sig_pep)
#
#        self.compare_files(actual_full_sequences_with_sigpep_fp,
#                           expected_full_sequences_with_sigpep_fp)
#
#        self.compare_files(actual_mature_seqs_fp, expected_mature_seq_fp)
#
#
#    def test_pipe_aborts_in_signalp_step_when_no_sigpeps(self):
#        """
#        Test execution abortion if there are no signal peptides detected
#        """
#
#        test_fas_no_sigpep = os.path.join(self.tmp_dir, 'test_no_sigpep.fas')
#        tmp_dir_to_test_removal = os.path.join(self.tmp_dir, 'intermediate_tmp')
#        os.mkdir(tmp_dir_to_test_removal)
#
#        with self.assertRaises(SystemExit):
#            _, _, _ = utils.signalp(test_fas_no_sigpep,
#                                    tmp_dir_to_test_removal,
#                                    '',
#                                    '')
#
#        self.assertFalse(os.path.exists(tmp_dir_to_test_removal))
#
#
#    def test_detect_and_output_transporters_via_signalp_when_no_sigpeps(self):
#        """
#        Test optional transporter prediction is generated even if
#        abortion of script due to no signal peptides
#        """
#
#        test_fas_no_sigpep = os.path.join(self.tmp_dir, 'test_no_sigpep.fas')
#        tmp_dir_to_test_removal = os.path.join(self.tmp_dir, 'intermediate_tmp_2')
#        os.mkdir(tmp_dir_to_test_removal)
#        rename_mappings = {'KGQ972831Candidaalbicans': 'KGQ972831Candidaalbicans'}
#
#        with self.assertRaises(SystemExit):
#            _, _, _ = utils.signalp(test_fas_no_sigpep,
#                                    tmp_dir_to_test_removal,
#                                    rename_mappings,
#                                    'actual',
#                                    trans=2)
#        self.assertFalse(os.path.exists(tmp_dir_to_test_removal))
#
#        expected_transporters = os.path.join(self.tmp_dir,
#                                             'expected_no_sigpep_transporters.fas')
#
#        actual_transporters = 'actual_predicted_transporters.fas'
#        self.assertTrue(actual_transporters)
#        self.compare_files(actual_transporters, expected_transporters)
#
#
#    def test_wolfpsort_func(self):
#        """
#        Test that wolfpsort's output is correctly parsed
#        """
#
#        actual_extracellular_accessions = utils.wolfpsort(self.test_fas)
#
#        self.assertEqual(["PM50_trimmed_paired_qual32_contig_10408",
#                          "R1qual32.paired_(paired)_contig_1567"],
#                         actual_extracellular_accessions)
#
#
#    def test_targetp_func(self):
#        """
#        Test that targetp output is correctly parsed by function
#        """
#
#        actual_secreted_accessions = utils.targetp(self.test_fas)
#
#        self.assertEqual(["PM50_trimmed_paired_",
#                          "R1qual32.paired_(pai",
#                          "PM50_trimmed_paired_",
#                          "PM50_trimmed_paired_",
#                          "PM30_NoIndex_L003_R1",
#                          "R1qual32.paired_(pai",
#                          "PM30_NoIndex_L003_R1"],
#                         actual_secreted_accessions)
#
#
#if __name__ == '__main__':
#    unittest.main(buffer=True)
