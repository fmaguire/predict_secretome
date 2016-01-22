import os
import warnings
import sys
import re
import shutil
from predict_secretome import utils
import hashlib
import subprocess
from Bio import SeqIO

class predictSecretome(object):
    """
    Runner class for predicting the secretome
    from a protein fasta input file
    """
    def __init__(self, argv):
        """
        Initialise the namespace mainly
        by:
        - converting user input to a self.settings dicts
        - initialising empty dict for predictions from each tools - self.outputs
        - initialise dict containing all the relevant file paths - self.files
        - generating directory structure for run
        """
        self.dep_path = os.path.join(os.getcwd(), 'predict_secretome', 'dependencies',
                                     'bin')

        self.settings = {'check_run': argv.check,
                         'force_overwrite': argv.force,
                         'permissive': argv.permissive,
                         'transporters': argv.transporters}


        # parse the typing argument
        if argv.plant:
            euk_type = 'plant'
        elif argv.animal:
            euk_type = 'animal'
        elif argv.fungi:
            euk_type = 'fungi'
        self.settings.update({'euk_type': euk_type})


        # run check dependencies util
        utils.check_dependencies(self.dep_path,
                                 self.settings['check_run'])

        self.files = {'input_fasta': argv.input_file}

        if not argv.run_name:
            basename = os.path.basename(self.files['input_fasta'])
            argv.run_name = os.path.splitext(basename)[0]

        self.run_name = argv.run_name

        if self.settings['permissive']:
            combination_flag = "permissive"
        else:
            combination_flag = "conservative"

        self.settings.update({'output_dir': os.path.join(os.getcwd(),
            '{0}_secretome_prediction_{1}_{2}_output'.format(self.run_name,
                                                             self.settings['euk_type'],
                                                             combination_flag))})

        self.settings.update({'work_dir': os.path.join(self.settings['output_dir'],
                                                           "intermediate_outputs")})

        self.outputs = {}

        self._initialise_directory_structure()



    def _initialise_directory_structure(self):
        """
        Initialise the output_structure checking
        if it already exists and whether to overwrite it.
        """
        path_exists = os.path.exists(self.settings['output_dir'])

        if path_exists and not self.settings['force_overwrite']:
            sys.exit('Output directory ({}) already exists.'
                ' Delete or re-run pipeline with --force option to overwrite'.format(self.settings['output_dir']))
        elif path_exists and self.settings['force_overwrite']:
            shutil.rmtree(self.settings['output_dir'])

        # makedirs recursively makes full path
        # i.e. outputdir and outputdir/workdir
        os.makedirs(self.settings['work_dir'])


    def _reformat_input(self):
        """
        Format input fasta and write to output:
        Specifically replace accessions with hash of name to guarantee under 20 chars
        input: input_file - filename
               output_file - filename
               path - path to dependencies bins
        output: rename_mappings - a dict of {new_acc: original_acc}
        """

        print("\n##Formatting input fasta: {0}##".format(self.files['input_fasta']))


        formatted_fasta_fp = os.path.join(self.settings['work_dir'],
                                          "formatted_input.fasta")

        # as some dependencies (tmhmm) don't support accessions over 20 chars
        # replace accessions with a 19-char truncated hash of original accession
        # keeping track in a dict where key is the hash and val is original accession
        raw_fasta_in_fh = open(self.files['input_fasta'], "r")
        formatted_out_fh = open(formatted_fasta_fp, 'w')

        # SeqIO doesn't expose the wrap attribute so need to call FastaIO directly
        # to use single line format
        fasta_out = SeqIO.FastaIO.FastaWriter(formatted_out_fh, wrap=None)
        fasta_out.write_header()
        rename_mappings = {}

        for record in SeqIO.parse(raw_fasta_in_fh, 'fasta'):
            truncated_md5 = hashlib.md5(\
                  record.description.encode('utf-8')).hexdigest()[:19]
                  # hash collision calculation fun time:
                  # md5 is 128bits (hexdigest 32 chars) so a 19 character truncation
                  # is 76 bits.  Therefore can calculate collision rate using birthday
                  # attack

            rename_mappings.update({truncated_md5: record.description})

            record.id = truncated_md5
            record.name = ''
            record.description = ''

            fasta_out.write_record(record)

        fasta_out.write_footer()

        raw_fasta_in_fh.close()
        formatted_out_fh.close()

        print("Input fasta formatted: {0}".format(formatted_fasta_fp))

        # update namespace on completion
        self.files.update({'reformatted_input': formatted_fasta_fp})
        self.rename_mappings = rename_mappings



    def _signalp(self):
        """
        Use signalp to identify sequences with signal peptides
        and create 2 files:
            - the full sequences of those with signal peptides
            - all input sequences with any signal peptides cleaved
        """

        # Dependency
        if 'reformatted_input' not in self.files.keys():
            self._reformat_input()


        # mature seqs are those with signal peptides removed
        mature_seqs_fp = os.path.join(self.settings['work_dir'],
                                      'signalp_mature_seqs.fasta')

        sigp = os.path.join(self.dep_path, 'signalp')
        sigp_cmd = "{0} -t euk -f short -m {1} {2}".format(sigp,
                                                           mature_seqs_fp,
                                                           self.files['reformatted_input'])


        print('##Detecing signal peptides##')
        sigp_proc = subprocess.Popen(sigp_cmd.split(),
                                     stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE)
        _, stderr = sigp_proc.communicate()
        print("Search Complete")

        # if there are no signal peps - then just use output from wolfpsort
        # and tmhmm filtering
        if '# No sequences predicted with a signal peptide' in stderr.decode('ascii'):
            print("No signal peptides detected so secretome can only be predicted using"
                   "permissive settings")
            return

        print("Collecting accessions with signal peptides")
        # get list of all accessions with signal peptides from mature
        # sequences file
        accessions_with_sig_pep = []
        with open(mature_seqs_fp, 'r') as mature_fh:
            for line in mature_fh.readlines():
                if line.startswith('>'):
                    strip_line = line.lstrip('>')
                    accessions_with_sig_pep.append(strip_line.split(' ; ')[0])


        # use this list to get full sequences (i.e. not mature) with signal peptides from
        # formatted input sequences
        full_sequences_with_sigpep_fp = os.path.join(self.settings['work_dir'],
                                                     'signalp_full_seqs_with_sigpep.fasta')

        print("Assembling full sequences with signal peptides")

        utils.write_seqs_from_accessions(accessions_with_sig_pep,
                                         self.files['reformatted_input'],
                                         full_sequences_with_sigpep_fp)
        print("Generating signal peptide free input fasta")
        # i.e. taking all sequences without signal peptides and the mature
        # sequences of those that do have them.

        # get sequences without signal peptides and write them to mature
        # seqs file

        accessions_without_sig_pep = [acc for acc in self.rename_mappings.keys()
                                        if acc not in accessions_with_sig_pep]

        input_without_sigpeps = os.path.join(self.settings['work_dir'],
                                             'input_seqs_without_sigpeps')

        # as we want to append the non signal peptide sequences
        # to the mature sequences
        shutil.copy(mature_seqs_fp, input_without_sigpeps)


        utils.write_seqs_from_accessions(accessions_without_sig_pep,
                                         self.files['reformatted_input'],
                                         input_without_sigpeps,
                                         append=True)


        self.files.update({'mature_sequences': mature_seqs_fp,
                           'sig_pep_free_seqs': input_without_sigpeps,
                           'sequences_with_sig': full_sequences_with_sigpep_fp})

        self.outputs.update({'acc_with_sig': accessions_with_sig_pep})


    def _tmhmm(self):
        """
        Run tmhmm on all sequences (if those sequences have signalpeptides
        then trim them off i.e. use mature seqs).

        Any sequences with TM domains are unlikely to be secreted.
        """

        #if "mature_seqs" not in self.files.keys():
        #    self._signalp()

        tmhmm = os.path.join(self.dep_path, 'tmhmm')
        tmhmm_cmd = "{0} {1}".format(tmhmm,
                                     self.files['sig_pep_free_seqs'])


        print("##Search for TM domains in sequences without ##")
        tmhmm_output = subprocess.check_output(tmhmm_cmd.split())
        tmhmm_output = tmhmm_output.decode('ascii').split('\n')
        print("Search complete")

        with open(os.path.join(self.settings['work_dir'], "tmhmm_output"), 'w') as tmhmm_fh:
            for line in tmhmm_output:
                tmhmm_fh.write(line + '\n')



        print("Parsing results")
        # Parse tmhmm raw output and write acc without tm domains
        # in non signal peptide sequence to output_file

        acc_without_tm_in_mature_seq = []
        putative_transporter = []

        for line in tmhmm_output:
            split_line = line.split('\t')
            if len(split_line) != 6:
                continue
            num_helices = int(split_line[4].split('=')[-1])
            if num_helices == 0:
                acc_without_tm_in_mature_seq.append(split_line[0])
            if self.settings['transporters'] and num_helices >= self.settings['transporters']:
                putative_transporter.append(split_line[0])

        self.outputs.update({'non_tm_mature_accs': acc_without_tm_in_mature_seq})

        if self.settings['transporters']:
            print("##Generating List of Putative Transporters##")
            transporters = os.path.join(self.settings['output_dir'],
                                        'transporters_{}_tm_domains'.format(self.settings['transporters']))

            utils.write_seqs_from_accessions(putative_transporter,
                                             self.files['reformatted_input'],
                                             transporters)

            self.files.update({'transporters': transporters})


    def _targetp(self):
        """
        Runs targetp on full seqs of those identified as having
        signal peptides by signalp to identify those designated
        as targeted 'secreted'
        input: mature_seqs_fp - fasta file without signal peptides
               plant - boolean to use plant settings or not for targetp
         /utput: secreted_accessions - list of accessions with signalpeps
                                      predicted to be secreted

        """
        if 'sequences_with_sig' not in self.files.keys():
            self._signalp()

        if self.settings['euk_type'] == 'plant':
            targetp_flag = "-P"
        else:
            targetp_flag = "-N"

        # the fortran 'How' ANN classifier dependency bundled with targetp
        # can't handle paths longer than 80 chars so as a hack fix I'm creating
        # a symlink to /home/user and removing it after execution
        username = os.getlogin()
        home_targetp = '/home/{0}/targetp-1.1'.format(username)

        try:
            os.symlink(os.path.join(os.sep.join(self.dep_path.split(os.sep)[:-1]),
                                    'targetp-1.1'),
                       home_targetp)
        except FileExistsError:
            print("Symlink hack for targetp already exists, proceeding anyway")
        # FileExistsError:
        targetp_path = os.path.join(self.dep_path, 'targetp')


        # so due to the problems with targetp there is some limit to the total
        # fasta input length therefore need to crop sequences and split fasta file
        print("Splitting and truncating signal peptide containing sequences for "
              "targetp")
        record_iterator = SeqIO.parse(open(self.files['sequences_with_sig']),
                                      'fasta')

        partial_fasta_with_sig_fps = []
        for i, batch in enumerate(utils.batch_iterator(record_iterator,
                                                       100)):
            split_filename = os.path.join(self.settings['work_dir'],
                                          "cropped_sig_partial_{}.fasta".format(i))

            with open(split_filename, 'w') as partial_fh:
                partial_fasta = SeqIO.FastaIO.FastaWriter(partial_fh, wrap=None)
                partial_fasta.write_header()

                for record in batch:
                    # truncate
                    if len(record.seq) > 500:
                        record.seq = record.seq[:500]
                    partial_fasta.write_record(record)
                partial_fasta.write_footer()

            partial_fasta_with_sig_fps.append(split_filename)

        # run targetp
        print("Running targetp on split files")
        targetp_output = []
        for subfile in partial_fasta_with_sig_fps:
            targetp_cmd = "{0} {1} {2}".format(targetp_path,
                                               targetp_flag,
                                               subfile)

            targetp_proc = subprocess.Popen(targetp_cmd.split(),
                                            stdout=subprocess.PIPE,
                                            stderr=subprocess.PIPE)
            targetp_stdout, targetp_stderr = targetp_proc.communicate()
            if 'Aborted' in targetp_stderr.decode('ascii'):
                print(targetp_cmd)
                print("targetp has crashed, unfortunately as the core of"
                       "targetp is only distributed as a trained and compiled fortran"
                       "ANN this is hard to debug.  Possibly your input file is too"
                       "large, attempt splitting and re-running using partial"
                       "files")
                sys.exit(1)

            partial_output = targetp_stdout.decode('ascii').split('\n')[8:-3]
            targetp_output = targetp_output + partial_output

        print("Targetp search complete")

        #os.unlink(home_targetp)


        with open(os.path.join(self.settings['work_dir'], "targetp_out"), 'w') as targetp_fh:
            for line in targetp_output:
                targetp_fh.write(line + '\n')


        print("Parsing results")
        # remove header and tail cruft in targetp output
        # and get those that have S as top predicted target
        secreted_accessions = []
        for line in targetp_output:
            split_line = re.split('\s+', line)
            if split_line[5] == 'S':
                secreted_accessions.append(split_line[0])
        self.outputs.update({'secreted_acc': secreted_accessions})

    def _wolfpsort(self):
        """
        Run wolfPsort on all formatted sequences using fungi setting
        to get a predicted list of 'extracellular' accessions
        input: formatted_fasta_fp - formatted fasta file of all input seqs
        output: extracellular_accessions - list of accessions with signalpeps
                                            predicted to be extracellular

        """

        if 'reformatted_input' not in self.files.keys():
            self._reformat_input()

        wolfpsort_path = os.path.join(self.dep_path, 'runWolfPsortSummary')

        wolfpsort_cmd = "{0} {1} < {2}".format(wolfpsort_path,
                                               self.settings['euk_type'],
                                               self.files['reformatted_input'])

        print("##Identifying sequences belonging to 'extracellular' compartment##")

        wolfpsort_output = subprocess.check_output(wolfpsort_cmd, shell=True)
        wolfpsort_output = wolfpsort_output.decode('ascii').split('\n')

        with open(os.path.join(self.settings['work_dir'], "wolfpsort_out"), 'w') as wolf_fh:
            for line in wolfpsort_output:
                wolf_fh.write(line + '\n')

        print("Search complete")


        print("Parsing results")
        # Removes header from output
        extracellular_accessions = []
        for line in wolfpsort_output[1:-1]:
            split_line = re.split('\s+', line)
            if split_line[1] == 'extr':
                extracellular_accessions.append(split_line[0])
        self.outputs.update({'extracellular_acc': extracellular_accessions})


    def _get_secretome_accessions(self):
        """
        Get secretome accessions from the predictions
        """
        # if predictions not done, then do them
        print("##Combining predictions##")

        sig_peptides = set(self.outputs['acc_with_sig'])
        no_tm_in_mature = set(self.outputs['non_tm_mature_accs'])
        secreted_sig = set(self.outputs['secreted_acc'])
        extracellular = set(self.outputs['extracellular_acc'])

        print("Accessions with signal peptides (signalp):", len(sig_peptides), sig_peptides)
        print("Accessions without a TM in mature sequence (tmhmm):", len(no_tm_in_mature), no_tm_in_mature)
        print("Secreted signal peptide (targetp):", len(secreted_sig), secreted_sig)
        print("Extracellular (wolfpsort):", len(extracellular), extracellular)

        # if conservative only get those accessions predicted as
        # A) having a signal peptide (signalp)
        # B) not having any TM domains outside of this sigpep
        #     (tmhmm)
        # C) a signal peptide predicted as being for secretion
        #     (targetp)
        # D) a sequence predicted as belonging to the extracellular
        #     compartment (wolFPSort)
        # if permissive get all accessions that belong to any of these
        # categories
        signal_pep_secreted = set.intersection(sig_peptides,
                                               no_tm_in_mature,
                                               secreted_sig)


        if not self.settings['permissive']:
            self.out_flag = 'conservative'
            predicted_acc_list = set.intersection(signal_pep_secreted,
                                                  extracellular)
        else:
            self.out_flag = 'permissive'
            # maybe remove no_tm_domains from this
            # as plenty of things don't have tm domains
            # that aren't secreted
            predicted_acc_list = set.union(signal_pep_secreted,
                                           extracellular)

        if len(predicted_acc_list) is 0:
            print("No secreted proteins found using {0} setting".format(self.out_flag))

        self.outputs.update({'predicted_secretome_acc': predicted_acc_list})



    def _output_secretome(self):
        """
        Using the secretome accessions retrieve sequences from
        the formatted renamed input fasta and reconstitute their
        names.

        Output file to output fasta.
        """

        print("##Writing predicted secretome fasta file##")

        formatted_input_fp = open(self.files['reformatted_input'], 'r')

        secretome_output = os.path.join(self.settings['output_dir'],
                                        "{0}_{1}_predicted_secretome.fasta".format(self.run_name,
                                                                                   self.out_flag))

        out_handle = open(secretome_output, 'w')

        fasta_out = SeqIO.FastaIO.FastaWriter(out_handle, wrap=None)
        fasta_out.write_header()

        records_written = 0

        print("Retreving Secretome fasta")

        for record in SeqIO.parse(formatted_input_fp, 'fasta'):
            if record.description in self.outputs['predicted_secretome_acc']:
                record.id = self.rename_mappings[record.id]
                record.name = ''
                record.description = ''
                fasta_out.write_record(record)
                records_written+=1

        if records_written > 0:
            fasta_out.write_footer()
        else:
            print("No secreted proteins found")
        out_handle.close()
        print("Secretome identification Complete")

    def run_predict(self):
        """
        Main runner method for the class that runs
        the series of steps required for generation
        of a predict secretome
        """
        self._reformat_input()
        self._signalp()
        self._tmhmm()
        self._targetp()
        self._wolfpsort()
        self._get_secretome_accessions()
        self._output_secretome()
        return
