import os
import pandas as pd
import copy as COPY
import matplotlib.pyplot as plt
import numpy as np
from .diachromatic_interaction_set import DiachromaticInteractionSet
from .ia_freq_dist_analysis import IaFreqDistAnalysis




class ReadpairAndInteractionCounter:
    """
    protocol = 'CHC' (capture Hi C) or 'HC' (Hi C)
    cell type = one of the abbreviations for one of the Javierre cell types, e.g., 'MK'
    min_interaction_distance: an integer such as 20000 for the minimum distance of an interaction to be included
    rpc_rule: e.g., 'ht'
    fdr: e.g., 0.05
    rpc_rule = 'ht'
    fdr = '05'

# Create DiachromaticInteractionSet
d11_interaction_set = DiachromaticInteractionSet(rpc_rule = RPC_RULE)
d11_interaction_set.parse_file(
    i_file = INTERACTION_FILE,
    verbose = True)

# Pass the DiachromaticInteractionSet to IaFreqDistAnalysis
ia_freq_dist_analysis = IaFreqDistAnalysis()
report_dict = ia_freq_dist_analysis.ingest_interaction_set(
    d11_inter_set = d11_interaction_set,
    verbose = True)
    """
    def __init__(self, cell_type, interaction_file, outprefix, min_interaction_dist=20000, protocol='CHC', rpc_rule="ht", fdr=0.05) -> None:
        allowable_cell_types = {'MAC','MK'}
        if not cell_type in allowable_cell_types:
            raise ValueError(f"Did not recognize Javier Cell Type {cell_type}")
        self._cell_type = cell_type
        self._min_interaction_dist = min_interaction_dist
        if not protocol in {'CHC', 'HC'}:
            raise ValueError(f"Invalid protocol {protocol}")
        self._protocol = protocol
        self._rpc_rule = rpc_rule
        self._fdr = fdr
        if not os.path.isfile(interaction_file):
            raise FileNotFoundError(f"Could not find file {interaction_file}")
        self._outprefix = outprefix
        
        d11_interaction_set = DiachromaticInteractionSet(rpc_rule = self._rpc_rule)
        d11_interaction_set.parse_file( i_file = interaction_file, verbose = True)
        self._ia_freq_dist_analysis = IaFreqDistAnalysis()
        self._report_dict = self._ia_freq_dist_analysis.ingest_interaction_set(d11_inter_set = d11_interaction_set,verbose = True)
        # set up dictionary for counts
        # we count the same items for each of the categories DIX, DI, UIR, UI, and ALL.
        my_count_dict = {
        'NN': {'S1': 0, 'S2': 0, 'T1': 0, 'T2': 0},
        'EE': {'S1': 0, 'S2': 0, 'T1': 0, 'T2': 0},
        'NE': {'S1': 0, 'S2': 0, 'T1': 0, 'T2': 0},
        'EN': {'S1': 0, 'S2': 0, 'T1': 0, 'T2': 0},
        'ALL': {'S1': 0, 'S2': 0, 'T1': 0, 'T2': 0}
        }
        RP_TYPE_FREQ_DICT = {
                'DIX': COPY.deepcopy(my_count_dict),
                 'DI': COPY.deepcopy(my_count_dict),
                'UIR': COPY.deepcopy(my_count_dict),
                 'UI': COPY.deepcopy(my_count_dict),
                'ALL': COPY.deepcopy(my_count_dict)
        }
        RP_TYPE_DENS_DICT = COPY.deepcopy(RP_TYPE_FREQ_DICT)
        for d11_inter in d11_interaction_set.interaction_list:
        # Get tags for interaction and enrichment category
            i_cat = d11_inter.get_category()
            e_cat = d11_inter.enrichment_status_tag_pair
    
            # Counts within each category
            RP_TYPE_FREQ_DICT[i_cat][e_cat]['S1'] += d11_inter._simple_1
            RP_TYPE_FREQ_DICT[i_cat][e_cat]['S2'] += d11_inter._simple_2
            RP_TYPE_FREQ_DICT[i_cat][e_cat]['T1'] += d11_inter._twisted_1
            RP_TYPE_FREQ_DICT[i_cat][e_cat]['T2'] += d11_inter._twisted_2
    
            # Counts combined for interaction categories
            RP_TYPE_FREQ_DICT['ALL'][e_cat]['S1'] += d11_inter._simple_1
            RP_TYPE_FREQ_DICT['ALL'][e_cat]['S2'] += d11_inter._simple_2
            RP_TYPE_FREQ_DICT['ALL'][e_cat]['T1'] += d11_inter._twisted_1
            RP_TYPE_FREQ_DICT['ALL'][e_cat]['T2'] += d11_inter._twisted_2
    
            # Counts combined for enrichment categories
            RP_TYPE_FREQ_DICT[i_cat]['ALL']['S1'] += d11_inter._simple_1
            RP_TYPE_FREQ_DICT[i_cat]['ALL']['S2'] += d11_inter._simple_2
            RP_TYPE_FREQ_DICT[i_cat]['ALL']['T1'] += d11_inter._twisted_1
            RP_TYPE_FREQ_DICT[i_cat]['ALL']['T2'] += d11_inter._twisted_2
    
            # Counts combined for interaction and enrichment categories
            RP_TYPE_FREQ_DICT['ALL']['ALL']['S1'] += d11_inter._simple_1
            RP_TYPE_FREQ_DICT['ALL']['ALL']['S2'] += d11_inter._simple_2
            RP_TYPE_FREQ_DICT['ALL']['ALL']['T1'] += d11_inter._twisted_1
            RP_TYPE_FREQ_DICT['ALL']['ALL']['T2'] += d11_inter._twisted_2

            # Fill second dictionary with relative frequencies
            
            for i_cat in ['DIX','DI','UIR','UI','ALL']:
                for e_cat in ['NN','EE','NE','EN','ALL']:
                    rp_total = sum(RP_TYPE_FREQ_DICT[i_cat][e_cat].values())
                    if 0 < rp_total:
                        RP_TYPE_DENS_DICT[i_cat][e_cat]['S1'] = RP_TYPE_FREQ_DICT[i_cat][e_cat]['S1']/rp_total
                        RP_TYPE_DENS_DICT[i_cat][e_cat]['S2'] = RP_TYPE_FREQ_DICT[i_cat][e_cat]['S2']/rp_total
                        RP_TYPE_DENS_DICT[i_cat][e_cat]['T1'] = RP_TYPE_FREQ_DICT[i_cat][e_cat]['T1']/rp_total
                        RP_TYPE_DENS_DICT[i_cat][e_cat]['T2'] = RP_TYPE_FREQ_DICT[i_cat][e_cat]['T2']/rp_total
                    else:
                        RP_TYPE_DENS_DICT[i_cat][e_cat]['S1'] = 0.0
                        RP_TYPE_DENS_DICT[i_cat][e_cat]['S2'] = 0.0
                        RP_TYPE_DENS_DICT[i_cat][e_cat]['T1'] = 0.0
                        RP_TYPE_DENS_DICT[i_cat][e_cat]['T2'] = 0.0
            self._RP_TYPE_FREQ_DICT = RP_TYPE_FREQ_DICT
            self._RP_TYPE_DENS_DICT = RP_TYPE_DENS_DICT
            
        my_ht_freq_dict = {
            'NN': {'0X': 0, '1X': 0,'2X': 0,'3X': 0, '01': 0,'02': 0,'03': 0,'12': 0,'13': 0,'23': 0},
            'EE': {'0X': 0, '1X': 0,'2X': 0,'3X': 0, '01': 0,'02': 0,'03': 0,'12': 0,'13': 0,'23': 0},
            'NE': {'0X': 0, '1X': 0,'2X': 0,'3X': 0, '01': 0,'02': 0,'03': 0,'12': 0,'13': 0,'23': 0},
            'EN': {'0X': 0, '1X': 0,'2X': 0,'3X': 0, '01': 0,'02': 0,'03': 0,'12': 0,'13': 0,'23': 0},
            'ALL': {'0X': 0, '1X': 0,'2X': 0,'3X': 0, '01': 0,'02': 0,'03': 0,'12': 0,'13': 0,'23': 0}
        }
        HT_TAG_FREQ_DICT = {
            'DIX': COPY.deepcopy(my_ht_freq_dict),
            'DI':  COPY.deepcopy(my_ht_freq_dict),
            'UIR': COPY.deepcopy(my_ht_freq_dict),
            'UI':  COPY.deepcopy(my_ht_freq_dict),
            'ALL': COPY.deepcopy(my_ht_freq_dict),
        }
        HT_TAG_DENS_DICT = COPY.deepcopy(HT_TAG_FREQ_DICT)
        # Get absolute frequencies
        for d11_inter in d11_interaction_set.interaction_list:
            i_cat = d11_inter.get_category()
            e_cat = d11_inter.enrichment_status_tag_pair
            ht_tag = d11_inter.get_ht_tag()
            HT_TAG_FREQ_DICT[i_cat][e_cat][ht_tag] += 1
            HT_TAG_FREQ_DICT['ALL'][e_cat][ht_tag] += 1
            HT_TAG_FREQ_DICT[i_cat]['ALL'][ht_tag] += 1
            HT_TAG_FREQ_DICT['ALL']['ALL'][ht_tag] += 1
    
        # Fill second dictionary with realtive frequencies
        for i_cat in ['DIX','DI','UIR','UI','ALL']:
            for e_cat in ['NN','EE','NE','EN','ALL']:
                i_total = sum(HT_TAG_FREQ_DICT[i_cat][e_cat].values())
                for ht_tag in ['0X','1X','2X','3X', '01','02','03','12','13','23']:
                    if 0 < i_total:
                        HT_TAG_DENS_DICT[i_cat][e_cat][ht_tag] = HT_TAG_FREQ_DICT[i_cat][e_cat][ht_tag]/i_total
                    else:
                        HT_TAG_DENS_DICT[i_cat][e_cat][ht_tag] = 0.0
        self._HT_TAG_FREQ_DICT = HT_TAG_FREQ_DICT
        self._HT_TAG_DENS_DICT = HT_TAG_DENS_DICT
        
        
    def get_interaction_set_summary(self):
        """_summary_
        

        return report
        """
        iset_d = self._ia_freq_dist_analysis._ingest_interaction_set_info_dict 
        entries = []
        for i_cat in ['DIX', 'DI', 'UIR', 'UI', 'ALL', ]:
            d = {'category': i_cat}
            total = 0
            for e_cat in ['NN', 'EE', 'NE', 'EN']:
                d[e_cat] = "{:,}".format(iset_d[i_cat][e_cat])
                total += iset_d[i_cat][e_cat]
            d['total'] = total
            entries.append(d)
        
        return pd.DataFrame(entries)
    
    def get_ht_tag_bar_chart_for_two_e_cats(self,
                                      e_cat_1 = 'NE',
                                      e_cat_2 = 'EN',
                                      e_cat_1_color = 'red',
                                      e_cat_2_color = 'blue',
                                        i_cats = ['DIX','DI','UIR','UI','ALL'],
                                      i_cat_colors = None,
                                      pdf_file_name = 'ht_tag_barplot_for_two_e_cats.pdf'
                                    ):
        """
        This function creates one grouped bar chart for each interaction category.
        For each HT tag there are two bars the represent two enrichment categories.
        """

        if i_cat_colors is None:
            i_cat_colors = {'DIX': 'orangered','DI': 'orange','UIR': 'lightblue','UI': 'lightgray','ALL': 'cornflowerblue'}

        # Fill second dictionary with recalcultated realtive frequencies for NE and EN combined (required for second y-axis)
        ht_tag_dens_dict = COPY.deepcopy(self._HT_TAG_FREQ_DICT)
        for i_cat in i_cats:
            i_total = sum(self._HT_TAG_FREQ_DICT[i_cat][e_cat_1].values()) + sum(self._HT_TAG_FREQ_DICT[i_cat][e_cat_2].values())
            
            for ht_tag in ['0X','1X','2X','3X', '01','02','03','12','13','23']:
                if 0 < i_total:
                    ht_tag_dens_dict[i_cat][e_cat_1][ht_tag] = self._HT_TAG_FREQ_DICT[i_cat][e_cat_1][ht_tag]/i_total
                    ht_tag_dens_dict[i_cat][e_cat_2][ht_tag] = self._HT_TAG_FREQ_DICT[i_cat][e_cat_2][ht_tag]/i_total
                else:
                    ht_tag_dens_dict[i_cat][e_cat_1][ht_tag] = 0.0
                    ht_tag_dens_dict[i_cat][e_cat_2][ht_tag] = 0.0

        # Determine maximal densities over all interactions categories and NE and EN
        y_max_d = 0.00
        for i_cat in i_cats:
            for e_cat in [e_cat_1, e_cat_2]:
                if y_max_d < max(ht_tag_dens_dict[i_cat][e_cat].values()):
                    y_max_d = max(ht_tag_dens_dict[i_cat][e_cat].values())
        y_max_d = y_max_d + y_max_d/10


        # Create a figure with plots for all interaction categories
        fig, ax = plt.subplots(5, figsize=(5,12))

        x_labels = ['0X','1X','2X','3X', '01','02','03','12','13','23']
        x = np.arange(len(x_labels))  # the label locations
        width = 0.35  # the width of the bars

        row_idx = 0 # Row in the plot grid
        for i_cat in i_cats:
            # The recalcultated densities for the combined enrichment categories
            ne_dens = ht_tag_dens_dict[i_cat][e_cat_1].values()
            en_dens = ht_tag_dens_dict[i_cat][e_cat_2].values()

            # The absolute frequencies that were passed to the function
            ne_freq = self._HT_TAG_FREQ_DICT[i_cat][e_cat_1].values()
            en_freq = self._HT_TAG_FREQ_DICT[i_cat][e_cat_2].values()

            # Create bars for densities
            rects_ne = ax[row_idx].bar(x - 0.5*width, ne_dens, width, label=e_cat_1, color=e_cat_1_color, edgecolor=i_cat_colors[i_cat], linewidth=3)
            rects_en = ax[row_idx].bar(x + 0.5*width, en_dens, width, label=e_cat_2, color=e_cat_2_color, edgecolor=i_cat_colors[i_cat], linewidth=3)

            # Add some text for labels, title and custom x-tick labels, etc.
            ax[row_idx].set_xlabel('HT tag')
            ax[row_idx].set_ylabel('Density')
            ax[row_idx].set_title(i_cat + ': Interaction numbers by HT tags', loc='left')
            ax[row_idx].set_xticks(x)
            ax[row_idx].set_xticklabels(x_labels)
            ax[row_idx].legend()

            # Add second axis with absolute frequencies
            ax2 = ax[row_idx].twinx()
            ax2.bar(x - 0.9*width, ne_freq, width, color='black', edgecolor='black', linewidth=3, alpha=0) # Set to alpha=1 to control scaling of y-axes
            ax2.bar(x + 0.9*width, en_freq, width, color='black', edgecolor='black', linewidth=3, alpha=0)
            ax2.set_ylabel('Frequency')

            # Scale second y-axis to maximum density
            ax[row_idx].set_ylim(0, y_max_d)    
            y_max_f = y_max_d*(sum(ne_freq) + sum(en_freq)) # This is possible, because we are using the combined densities of NE and EN
            ax2.set_ylim(0, y_max_f)
            #ax2.ticklabel_format(axis='y', style='sci', scilimits=(0,1))

            # Go to the next row of the plot grid
            row_idx += 1

        # Save and show plot
        fig.tight_layout()
        fig.savefig(pdf_file_name)
        plt.show()
    
    
    def get_readcount_freq_df(self):
        e_cat = 'ALL'
        entries = []
        for i_cat in ['DIX','DI','UIR','UI','ALL']:
            title = f"{i_cat}-{e_cat}"
            freq_d = self._RP_TYPE_FREQ_DICT[i_cat][e_cat]
            dens_d = self._RP_TYPE_DENS_DICT[i_cat][e_cat]
            d = {'title':title, "simple 1 (n)": str(freq_d['S1']), "simple 2 (n)": str(freq_d['S2']), 
                                "twisted 1 (n)": str(freq_d['T1']),"twisted 2 (n)": str(freq_d['T2']),
                                "simple 1 (%%)": "{:.2f}%".format(100*dens_d['S1']), 
                                "simple 2 (%%)": "{:.2f}%".format(100*dens_d['S2']),
                                "twisted 1 (%%)": "{:.2f}%".format(100*dens_d['T1']),
                                "twisted 2 (%%)": "{:.2f}%".format(100*dens_d['T2'])}
            entries.append(d)
            
        return pd.DataFrame(entries)
    
    
    def plot_category_distribution(self, pdf_file_name=None):
        e_cat = 'ALL'

        RP_CAT_COLORS = [(255/255, 160/255, 200/255), (255/255, 80/255, 120/255), (80/255, 190/255, 120/255), (60/255, 150/255, 120/255)]
        I_CAT_COLORS = {
            'DIX': 'orangered',
            'DI': 'orange',
            'UIR': 'lightblue',
            'UI': 'lightgray',
            'ALL': 'cornflowerblue'
        }
        rp_cat_labels = ['0', '1', '2', '3']

        # Determine y_max
        y_max_d = 0.25
        for i_cat in ['DIX', 'DI', 'UIR', 'UI', 'ALL']:
            if y_max_d < max(self._RP_TYPE_DENS_DICT[i_cat][e_cat].values()):
                y_max_d = max(self._RP_TYPE_DENS_DICT[i_cat][e_cat].values())
        y_max_d = 1.1* y_max_d  # add 10% extra space


        x = np.arange(len(rp_cat_labels))  # the label locations
        width = 0.35  # the width of the bars
        density_tick_labels = np.arange(0, 1, 0.25)

        fig, ax = plt.subplots(5, figsize=(5,12))

        row_idx = 0
        for i_cat in ['DIX', 'DI', 'UIR', 'UI', 'ALL']:
            # Create barchart
            ax[row_idx].bar(x, self._RP_TYPE_DENS_DICT[i_cat][e_cat].values(), width, color=RP_CAT_COLORS)
            ax[row_idx].set_title(i_cat + '-' + e_cat, loc='left')
            ax[row_idx].set_xticks(x)
            ax[row_idx].set_xticklabels(rp_cat_labels)
            ax[row_idx].axhline(0.25, zorder=0, color='gray', linewidth=0.5)
            ax[row_idx].set_xlabel('Read pair type')
            ax[row_idx].set_ylabel('Density')
            ax[row_idx].set_yticks(ticks=density_tick_labels)
    
            # Add second axis with absolute counts and normalize to maximum density
            ax2 = ax[row_idx].twinx()
            ax2.bar(x, self._RP_TYPE_FREQ_DICT[i_cat][e_cat].values(), width, color=RP_CAT_COLORS, edgecolor=I_CAT_COLORS[i_cat], linewidth=3)
            ax2.set_ylabel('Frequency')
            ax[row_idx].set_ylim(0, y_max_d)    
            y_lim = y_max_d*sum(self._RP_TYPE_FREQ_DICT[i_cat][e_cat].values())
            ax2.set_ylim(0, y_lim)

            row_idx += 1

            fig.tight_layout()
            if not pdf_file_name is None:
                fig.savefig(pdf_file_name)

        plt.show()

    def get_ht_summary_df(self):
        entries = []
        for i_cat in ['DIX', 'DI', 'UIR', 'UI', 'ALL']:
            for ht_tag in ['0X','1X','2X','3X', '01','02','03','12','13','23']:
                d = {'category' : i_cat, 'ht_tag': ht_tag}
                for e_cat in ['NN', 'EE', 'NE', 'EN', 'ALL']:
                    d[e_cat] = "{:,}".format(self._HT_TAG_FREQ_DICT[i_cat][e_cat][ht_tag])
                    perc = f"{e_cat} (%)"  
                    d[perc] =  "{:.2f}%".format(100*self._HT_TAG_DENS_DICT[i_cat][e_cat][ht_tag])                      
                entries.append(d) 
            
        return pd.DataFrame(entries)
        
        