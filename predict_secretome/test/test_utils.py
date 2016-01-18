#!/usr/bin/env python3

from predict_secretome.test.base import CoreTestingClass
from predict_secretome import utils

import os

class TestUtils(CoreTestingClass):
    """
    Class for testing the general utility functions
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



class TestFormatFasta(CoreTestingClass):
    """
    Class to test the fasta reformatting works and doesn't destroy information
    """

    def setUp(self):
        """
        Setup tests with shared input variables
        """
        self.test_fas = os.path.join(self.test_dir,
                                     "test.fas")

        self.expected_formatted_fasta = os.path.join(self.test_dir,
                                                     'expected_formatted_fasta')

    def test_format_fasta(self):
        """
        Test that the input fasta is correctly reformatted
        """

        tmp_dir = self.test_dir

        mappings, formatted_fasta = utils.format_fasta(self.test_fas, tmp_dir)

        self.compare_files(formatted_fasta, self.expected_formatted_fasta)

        self.assertEqual(len(mappings), 10)


    def test_renaming_fasta(self):
        """
        Test that output fasta is correctly renamed to original accessions
        """
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


