from unittest import TestCase
import os
import sys
from numpy import arange, log
import warnings
from diachr.diachromatic_interaction_set import DiachromaticInteractionSet
from diachr.randomize_interaction_set import RandomizeInteractionSet

PACKAGE_PARENT = '..'
SCRIPT_DIR = os.path.dirname(os.path.realpath(os.path.join(os.getcwd(), os.path.expanduser(__file__))))
sys.path.append(os.path.normpath(os.path.join(SCRIPT_DIR, PACKAGE_PARENT)))


class TestRandomizationAnalysis(TestCase):

    @classmethod
    def setUpClass(cls):

        # Mute warnings
        warnings.simplefilter('ignore')
        test_dir = os.path.dirname(__file__)

        # Prepare test data to test for correct use of P-value thresholds
        cls.interaction_set_1 = DiachromaticInteractionSet(rpc_rule='st')
        fdr_1_file = os.path.join(test_dir, "data/test_03/diachromatic_interaction_file_fdr_1.tsv.gz")
        cls.interaction_set_1.parse_file(fdr_1_file)
        cls.randomize_1 = RandomizeInteractionSet(random_seed=42)

        # Prepare interaction set using the test file with 1,000 interactions
        cls.interaction_set_1000 = DiachromaticInteractionSet(rpc_rule='st')
        fdr_top_1000_file = os.path.join(test_dir, "data/test_03/diachromatic_interaction_file_fdr_top_1000.tsv.gz")
        cls.interaction_set_1000.parse_file(fdr_top_1000_file)
        cls.randomize_1000 = RandomizeInteractionSet(random_seed=0)

        # Prepare interaction set using the test file with 64,000 interactions
        cls.interaction_set_64000 = DiachromaticInteractionSet(rpc_rule='st')
        fdr_top_64000_file = os.path.join(test_dir, "data/test_03/diachromatic_interaction_file_fdr_top_64000.tsv.gz")
        cls.interaction_set_64000.parse_file(fdr_top_64000_file)
        cls.randomize_64000 = RandomizeInteractionSet(random_seed=0)

        # Prepare interaction set to test for correct determination of potentially significant interaction numbers
        cls.interaction_set_pot_sig = DiachromaticInteractionSet(rpc_rule='st')
        fdr_pot_sig_file = os.path.join(test_dir, "data/test_03/diachromatic_interaction_file_test_pot_sig.tsv.gz")
        cls.interaction_set_pot_sig.parse_file(fdr_pot_sig_file)
        cls.randomize_pot_sig = RandomizeInteractionSet(random_seed=0)

    def test_parallel_processing(self):
        """
        Here it is tested whether the results are independent of the degree of parallel processing.
        """

        # Parameters for randomization
        nominal_alpha = 0.0025
        iter_num = 10

        # Perform randomization analysis without 'multiprocessing' package
        random_analysis_thread_num_0_info_dict = self.randomize_64000.perform_randomization_analysis(
            interaction_set=self.interaction_set_64000,
            nominal_alphas=[nominal_alpha],
            iter_num=iter_num,
            thread_num=0)

        # Perform randomization analysis with 'multiprocessing' package but only in one process
        random_analysis_thread_num_1_info_dict = self.randomize_64000.perform_randomization_analysis(
            interaction_set=self.interaction_set_64000,
            nominal_alphas=[nominal_alpha],
            iter_num=iter_num,
            thread_num=1)

        # Perform randomization analysis with 'multiprocessing' package in two parallel processes
        random_analysis_thread_num_2_info_dict = self.randomize_64000.perform_randomization_analysis(
            interaction_set=self.interaction_set_64000,
            nominal_alphas=[nominal_alpha],
            iter_num=iter_num,
            thread_num=2)

        # Perform randomization analysis with 'multiprocessing' package in three parallel processes
        random_analysis_thread_num_3_info_dict = self.randomize_64000.perform_randomization_analysis(
            interaction_set=self.interaction_set_64000,
            nominal_alphas=[nominal_alpha],
            iter_num=iter_num,
            thread_num=3)

        # Compare results
        for iter_idx in range(0, iter_num):
            sig_num_r_0 = random_analysis_thread_num_0_info_dict['SIG_NUM_R_LISTS'][nominal_alpha][iter_idx]
            sig_num_r_1 = random_analysis_thread_num_1_info_dict['SIG_NUM_R_LISTS'][nominal_alpha][iter_idx]
            sig_num_r_2 = random_analysis_thread_num_2_info_dict['SIG_NUM_R_LISTS'][nominal_alpha][iter_idx]
            sig_num_r_3 = random_analysis_thread_num_3_info_dict['SIG_NUM_R_LISTS'][nominal_alpha][iter_idx]
            self.assertEqual(sig_num_r_0, sig_num_r_1,
                             msg='Different numbers of randomized significant interactions depending on '
                                 'parallel processing!')
            self.assertEqual(sig_num_r_0, sig_num_r_2,
                             msg='Different numbers of randomized significant interactions depending on '
                                 'parallel processing!')
            self.assertEqual(sig_num_r_0, sig_num_r_3,
                             msg='Different numbers of randomized significant interactions depending on '
                                 'parallel processing!')

    def test_for_changes_in_results(self):
        """
        Here it is tested whether the results for a certain input have remain unchanged.
        """

        random_analysis_info_dict = self.randomize_64000.perform_randomization_analysis(
            interaction_set=self.interaction_set_64000,
            nominal_alphas=[0.0025],
            iter_num=100,
            thread_num=0)

        # Find index of the passed nominal alpha
        nominal_alpha_idx = random_analysis_info_dict['RESULTS']['NOMINAL_ALPHA'].index(0.0025)

        # The following results were obtained with an earlier version
        expected_sig_num_r_mean = 152.9
        expected_sig_num_r_sd = 12.47998
        expected_z_score = 101.85109

        # The following results are obtained with the current version
        observed_sig_num_r_mean = random_analysis_info_dict['RESULTS']['SIG_NUM_R_MEAN'][nominal_alpha_idx]
        observed_sig_num_r_sd = random_analysis_info_dict['RESULTS']['SIG_NUM_R_SD'][nominal_alpha_idx]
        observed_z_score = random_analysis_info_dict['RESULTS']['Z_SCORE'][nominal_alpha_idx]

        # Compare previous and current results
        self.assertAlmostEqual(expected_sig_num_r_mean, observed_sig_num_r_mean, places=5,
                               msg='A different mean number of randomized significant interactions was obtained for an '
                                   'earlier version!')
        self.assertAlmostEqual(expected_sig_num_r_sd, observed_sig_num_r_sd, places=5,
                               msg='A different standard deviation for the number of randomized significant '
                                   'interactions was '
                                   'obtained for an earlier version!')
        self.assertAlmostEqual(expected_z_score, observed_z_score, places=5,
                               msg='A different Z-score was obtained for an earlier version!')

    def test_results_for_same_and_different_random_seeds(self):
        """
        Here it is tested whether the same results are obtained with the same random seed.and whether different results
        are obtained with different random seeds.
        """

        # Perform analysis for RandomizeInteractionSet object of this class (random seed: 0)
        random_analysis_info_dict = self.randomize_1000.perform_randomization_analysis(
            interaction_set=self.interaction_set_1000,
            nominal_alphas=[0.005],
            iter_num=100,
            thread_num=0)

        # Create a new RandomizeInteractionSet object with the same random seed and perform the same analysis
        randomize_1000_same_seed = RandomizeInteractionSet(random_seed=0)
        random_analysis_info_dict_same_seed = randomize_1000_same_seed.perform_randomization_analysis(
            interaction_set=self.interaction_set_1000,
            nominal_alphas=[0.005],
            iter_num=100,
            thread_num=0)

        # Create a new RandomizeInteractionSet object with the different random seed and perform the same analysis
        randomize_1000_different_seed = RandomizeInteractionSet(random_seed=42)
        random_analysis_info_dict_different_seed = randomize_1000_different_seed.perform_randomization_analysis(
            interaction_set=self.interaction_set_1000,
            nominal_alphas=[0.005],
            iter_num=100,
            thread_num=0)

        # Find indices of the passed nominal alphas
        nominal_alpha_idx = random_analysis_info_dict['RESULTS']['NOMINAL_ALPHA'].index(0.005)
        nominal_alpha_idx_same_seed = random_analysis_info_dict_same_seed['RESULTS']['NOMINAL_ALPHA'].index(
            0.005)
        nominal_alpha_idx_different_seed = random_analysis_info_dict_different_seed['RESULTS'][
            'NOMINAL_ALPHA'].index(0.005)

        # Compare mean numbers of significant randomized interactions from the three objects
        sig_num_r_mean = random_analysis_info_dict['RESULTS']['SIG_NUM_R_MEAN'][nominal_alpha_idx]
        sig_num_r_mean_same_seed = random_analysis_info_dict_same_seed['RESULTS']['SIG_NUM_R_MEAN'][
            nominal_alpha_idx_same_seed]
        sig_num_r_mean_different_seed = \
            random_analysis_info_dict_different_seed['RESULTS']['SIG_NUM_R_MEAN'][nominal_alpha_idx_different_seed]
        self.assertEqual(sig_num_r_mean, sig_num_r_mean_same_seed, msg='Get different results with the same random '
                                                                       'seed!')
        self.assertNotEqual(sig_num_r_mean, sig_num_r_mean_different_seed, msg='Get the same results with different '
                                                                               'random seeds!')

    def test_independence_of_randomization_and_parallelization(self):
        """
        Here it is tested whether the same results are obtained for a given random seed, regardless of how many
        processes are used.
        """

        # Perform analysis without multiprocessing (random seed: 0)
        random_analysis_info_dict_0 = self.randomize_1000.perform_randomization_analysis(
            interaction_set=self.interaction_set_1000,
            nominal_alphas=[0.005],
            iter_num=100,
            thread_num=0)

        # Perform analysis with multiprocessing but only one process (random seed: 0)
        random_analysis_info_dict_1 = self.randomize_1000.perform_randomization_analysis(
            interaction_set=self.interaction_set_1000,
            nominal_alphas=[0.005],
            iter_num=100,
            thread_num=1)

        # Perform analysis with multiprocessing but only one process (random seed: 0)
        random_analysis_info_dict_2 = self.randomize_1000.perform_randomization_analysis(
            interaction_set=self.interaction_set_1000,
            nominal_alphas=[0.005],
            iter_num=100,
            thread_num=2)

        # Find indices of the passed nominal alphas
        nominal_alpha_idx_0 = random_analysis_info_dict_0['RESULTS']['NOMINAL_ALPHA'].index(0.005)
        nominal_alpha_idx_1 = random_analysis_info_dict_1['RESULTS']['NOMINAL_ALPHA'].index(0.005)
        nominal_alpha_idx_2 = random_analysis_info_dict_2['RESULTS']['NOMINAL_ALPHA'].index(0.005)

        # Compare mean numbers of significant randomized interactions from the three randomizations
        sig_num_r_mean_0 = random_analysis_info_dict_0['RESULTS']['SIG_NUM_R_MEAN'][nominal_alpha_idx_0]
        sig_num_r_mean_1 = random_analysis_info_dict_1['RESULTS']['SIG_NUM_R_MEAN'][nominal_alpha_idx_1]
        sig_num_r_mean_2 = random_analysis_info_dict_2['RESULTS']['SIG_NUM_R_MEAN'][nominal_alpha_idx_2]
        self.assertEqual(sig_num_r_mean_0, sig_num_r_mean_1, msg='Get different results with the same random depending '
                                                                 'multiprocessing!')
        self.assertEqual(sig_num_r_mean_0, sig_num_r_mean_2, msg='Get different results with the same random depending '
                                                                 'multiprocessing!')

    def test_determine_significant_pvals_at_nominal_alphas_nnl(self):
        """
        Here the function that determines the number of significant P-values for a list of nominal alphas is tested on
        a small example.
        """

        # Create list of nominal alphas
        nominal_alphas = [0.01, 0.02, 0.03, 0.04, 0.05]
        nominal_alphas_nnl = -log(nominal_alphas)

        # Create list of P-values
        p_values = [0.01, 0.011, 0.02, 0.02, 0.03, 0.03, 0.04, 0.041, 0.05, 0.05]
        p_values_nnl = -log(p_values)

        # Test the function implemented inn class RandomizeInteractionSet
        randomize = RandomizeInteractionSet()
        sig_num_list_nnl = randomize._determine_significant_pvals_at_nominal_alphas_log10(
            log10_nominal_alphas=nominal_alphas_nnl,
            log10_p_values=p_values_nnl)

        # Compare the returned list with the expected list
        expected_list = [1, 4, 6, 7, 10]
        for idx in range(0, len(sig_num_list_nnl)):
            self.assertEqual(expected_list[idx], sig_num_list_nnl[idx])

    def test_results_depending_on_how_a_nominal_alpha_was_passed(self):
        """
        Here it is tested whether the results are independent of the lists in which a certain nominal alpha is passed.

        Note: The largest nominal alpha determines the size of 'RP_INTER_DICT' and thus the randomization. In order to
        keep the randomization constant, the largest nominal alpha must be the same in all lists.
        """

        # Create lists of nominal alphas
        nominal_alphas_0 = [0.01, 0.02, 0.05]
        nominal_alphas_1 = [0.01, 0.05]
        nominal_alphas_2 = [0.02, 0.05]
        nominal_alphas_3 = [0.05]

        # Pass all nominal alphas together
        random_analysis_info_dict_0 = self.randomize_64000.perform_randomization_analysis(
            interaction_set=self.interaction_set_64000,
            nominal_alphas=nominal_alphas_0,
            iter_num=10,
            thread_num=0)

        # Pass only the first and last nominal alpha of the original list
        random_analysis_info_dict_1 = self.randomize_64000.perform_randomization_analysis(
            interaction_set=self.interaction_set_64000,
            nominal_alphas=nominal_alphas_1,
            iter_num=10,
            thread_num=0)

        # Pass only the second and last nominal alpha of the original list
        random_analysis_info_dict_2 = self.randomize_64000.perform_randomization_analysis(
            interaction_set=self.interaction_set_64000,
            nominal_alphas=nominal_alphas_2,
            iter_num=10,
            thread_num=0)

        # Pass only the third nominal alpha of the original list
        random_analysis_info_dict_3 = self.randomize_64000.perform_randomization_analysis(
            interaction_set=self.interaction_set_64000,
            nominal_alphas=nominal_alphas_3,
            iter_num=10,
            thread_num=0)

        # Compare Z-scores for first nominal alpha of the original list
        nominal_alpha_idx = random_analysis_info_dict_0['RESULTS']['NOMINAL_ALPHA'].index(nominal_alphas_0[0])
        z_score_0_0 = random_analysis_info_dict_0['RESULTS']['Z_SCORE'][nominal_alpha_idx]
        nominal_alpha_idx = random_analysis_info_dict_1['RESULTS']['NOMINAL_ALPHA'].index(nominal_alphas_0[0])
        z_score_1_0 = random_analysis_info_dict_1['RESULTS']['Z_SCORE'][nominal_alpha_idx]
        self.assertEqual(z_score_0_0, z_score_1_0, msg='Results differ depending on how a nominal alpha was passed!')

        # Compare Z-scores for second nominal alpha of the original list
        nominal_alpha_idx = random_analysis_info_dict_0['RESULTS']['NOMINAL_ALPHA'].index(nominal_alphas_0[1])
        z_score_0_1 = random_analysis_info_dict_0['RESULTS']['Z_SCORE'][nominal_alpha_idx]
        nominal_alpha_idx = random_analysis_info_dict_2['RESULTS']['NOMINAL_ALPHA'].index(nominal_alphas_0[1])
        z_score_2_1 = random_analysis_info_dict_2['RESULTS']['Z_SCORE'][nominal_alpha_idx]
        self.assertEqual(z_score_0_1, z_score_2_1, msg='Results differ depending on how a nominal alpha was passed!')

        # Compare Z-scores for third nominal alpha of the original list
        nominal_alpha_idx = random_analysis_info_dict_0['RESULTS']['NOMINAL_ALPHA'].index(nominal_alphas_0[2])
        z_score_0_2 = random_analysis_info_dict_0['RESULTS']['Z_SCORE'][nominal_alpha_idx]
        nominal_alpha_idx = random_analysis_info_dict_1['RESULTS']['NOMINAL_ALPHA'].index(nominal_alphas_0[2])
        z_score_1_2 = random_analysis_info_dict_1['RESULTS']['Z_SCORE'][nominal_alpha_idx]
        nominal_alpha_idx = random_analysis_info_dict_2['RESULTS']['NOMINAL_ALPHA'].index(nominal_alphas_0[2])
        z_score_2_2 = random_analysis_info_dict_2['RESULTS']['Z_SCORE'][nominal_alpha_idx]
        nominal_alpha_idx = random_analysis_info_dict_3['RESULTS']['NOMINAL_ALPHA'].index(nominal_alphas_0[2])
        z_score_3_2 = random_analysis_info_dict_3['RESULTS']['Z_SCORE'][nominal_alpha_idx]
        self.assertEqual(z_score_0_2, z_score_1_2, msg='Results differ depending on how a nominal alpha was passed!')
        self.assertEqual(z_score_0_2, z_score_2_2, msg='Results differ depending on how a nominal alpha was passed!')
        self.assertEqual(z_score_0_2, z_score_3_2, msg='Results differ depending on how a nominal alpha was passed!')

    def test_interactions_numbers_at_different_pval_thresholds(self):
        """
        Here it is tested whether the P-value thresholds are used correctly during the FDR procedure.
        """

        # Prepare list of nominal alphas
        nominal_alpha_max = 0.05
        nominal_alpha_step = 0.00025
        nominal_alphas = arange(nominal_alpha_step, nominal_alpha_max + nominal_alpha_step, nominal_alpha_step)

        # The randomization function returns a dictionary containing the numbers that are being tested here
        randomization_info_dict = self.randomize_1.perform_randomization_analysis(
            interaction_set=self.interaction_set_1,
            nominal_alphas=nominal_alphas,
            iter_num=1,
            thread_num=0,
            verbose=False)

        # We test the number of significant interactions for all P-value thresholds
        sig_num_o_expected = 10
        for sig_num_o in randomization_info_dict['RESULTS']['SIG_NUM_O']:
            self.assertEqual(sig_num_o_expected, sig_num_o,
                             msg='Number of significant interactions did not increase by 10!')
            sig_num_o_expected += 10

    def test_for_changes_in_fdr_results_top_64000_interactions(self):
        """
        Here it is tested whether the results of the randomization procedure for the test file with 64,000 input
        interactions have changed with respect to the determined P-value threshold that yields an FDR below 0.05.
        """

        # Prepare list of nominal alphas
        nominal_alpha_max = 0.05
        nominal_alpha_step = 0.00025
        nominal_alphas = arange(nominal_alpha_step, nominal_alpha_max + nominal_alpha_step, nominal_alpha_step)

        # The randomization function returns a dictionary containing the numbers that are being tested here
        randomization_info_dict = self.randomize_64000.perform_randomization_analysis(
            interaction_set=self.interaction_set_64000,
            nominal_alphas=nominal_alphas,
            iter_num=1,
            thread_num=0,
            verbose=False)

        # Index of the determined P-value threshold
        result_index = self.randomize_64000.get_largest_nominal_alpha_index_at_chosen_fdr_thresh(0.05)

        # The determined P-value threshold must be 0.00375
        determined_pval_thresh = randomization_info_dict['RESULTS']['NOMINAL_ALPHA'][result_index]
        self.assertAlmostEqual(0.00075, determined_pval_thresh, 5,
                               msg='When this test was created, a different P-value threshold was determined!')

        # At the determined P-value threshold the number of significant interactions must be 1,585
        determined_sig_num_o = randomization_info_dict['RESULTS']['SIG_NUM_O'][result_index]
        self.assertEqual(1004, determined_sig_num_o,
                         msg='When this test was created, a different number of significant interactions was '
                             'determined!')

        # At the determined P-value threshold the number of randomized significant interactions must be 77
        determined_sig_num_r = randomization_info_dict['RESULTS']['SIG_NUM_R_MEAN'][result_index]
        self.assertEqual(37, determined_sig_num_r,
                         msg='When this test was created, a different number of randomized significant interactions '
                             'was determined!')

        # At the determined P-value threshold the FDR must be 0.04858
        determined_fdr = randomization_info_dict['RESULTS']['FDR'][result_index]
        self.assertAlmostEqual(0.03685, determined_fdr, 5, msg='When this test was created, a different FDR was '
                                                               'determined!')

        # At the determined P-value threshold the Z-Score must be "NA"
        determined_z_score = randomization_info_dict['RESULTS']['Z_SCORE'][result_index]
        self.assertEqual("NA", determined_z_score, msg='When this test was created, a different Z-score was '
                                                       'determined!')

    def test_for_correct_determination_potentially_significant_interactions_numbers(self):
        """
        Here it is tested whether the number of potentially significant interactions at different nominal alphas is
        determined correctly. We prepared a test file that contains a total of nine interactions, three each with 7, 6
        and 5 read pairs. For each read pair number, there is one interaction with zero simple read pairs.
        """

        nominal_alpha_max = 0.07 # The smallest P-value with 5 read pairs (0.0625) is included
        nominal_alpha_step = 0.00001
        nominal_alphas = arange(nominal_alpha_step, nominal_alpha_max + nominal_alpha_step, nominal_alpha_step)

        # Perform randomization with the prepared interaction set
        randomize_pot_sig_info_dict = self.randomize_pot_sig.perform_randomization_analysis(
            interaction_set = self.interaction_set_pot_sig,
            nominal_alphas = nominal_alphas,
            iter_num=1)

        # Check interaction numbers at the transitions to largerr nominal alphas and fewer required read pairs
        nominal_alphas = randomize_pot_sig_info_dict['RESULTS']['NOMINAL_ALPHA']

        # Up to a nominal alpha of 0.03125, at least 7 read pairs are required for significance
        nominal_alpha_idx = nominal_alphas.index(0.03125)
        # We have three interactions with 7 read pairs
        pot_sig_num = randomize_pot_sig_info_dict['RESULTS']['POT_SIG_NUM'][nominal_alpha_idx]
        self.assertEqual(3, pot_sig_num)
        # One of these interactions has zero simple read pairs and is significant at a nominal alpha of 0.03125
        sig_num_o = randomize_pot_sig_info_dict['RESULTS']['SIG_NUM_O'][nominal_alpha_idx]
        self.assertEqual(1, sig_num_o)
        # For the next larger nominal alpha, only 6 read pairs are required
        nominal_alpha_idx = nominal_alphas.index(0.03126)
        # and we have six potentially significant interactions
        pot_sig_num = randomize_pot_sig_info_dict['RESULTS']['POT_SIG_NUM'][nominal_alpha_idx]
        self.assertEqual(6, pot_sig_num)
        # and two significant interactions
        sig_num_o = randomize_pot_sig_info_dict['RESULTS']['SIG_NUM_O'][nominal_alpha_idx]
        self.assertEqual(2, sig_num_o)

        # Up to a nominal alpha of 0.06250, at least 6 read pairs are required for significance
        nominal_alpha_idx = nominal_alphas.index(0.06250)
        # We have six interactions with at least 6 read pairs
        pot_sig_num = randomize_pot_sig_info_dict['RESULTS']['POT_SIG_NUM'][nominal_alpha_idx]
        self.assertEqual(6, pot_sig_num)
        # Two of these interactions have zero simple read pairs and are significant at a nominal alpha of 0.06250
        sig_num_o = randomize_pot_sig_info_dict['RESULTS']['SIG_NUM_O'][nominal_alpha_idx]
        self.assertEqual(2, sig_num_o)
        # For the next larger nominal alpha, only 5 read pairs are required
        nominal_alpha_idx = nominal_alphas.index(0.06251)
        # and we have nine potentially significant interactions
        pot_sig_num = randomize_pot_sig_info_dict['RESULTS']['POT_SIG_NUM'][nominal_alpha_idx]
        self.assertEqual(9, pot_sig_num)
        # and three significant interactions
        sig_num_o = randomize_pot_sig_info_dict['RESULTS']['SIG_NUM_O'][nominal_alpha_idx]
        self.assertEqual(3, sig_num_o)
