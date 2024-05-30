# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 09:39:39 2023

@author: Administrator
"""
import numpy as np
import matplotlib.pyplot as plt
from ric_fun import extract_matrix
from plots import joint_heatmap_hist, joint_heatmap_hist_genome
from read_genome_sites import dict_gene_sites_sars1, dict_gene_sites_ibv, dict_gene_sites_pdcov, dict_gene_sites_pedv, dict_gene_sites_sars, dict_gene_sites_mers
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns

res = 50
dir_pre = 'G:/'
file_dir = dir_pre + 'OneDrive - webmail.hzau.edu.cn/vRic-seq/data/rna-seq/' #'C:/Users/61918/Desktop/'
# ma_pedv_rna_cell_1 = extract_matrix(file_dir, 'PEDV_rna-cell-1.1nt.none.matrix', res)[0] 
ma_pedv_rna_cell = extract_matrix(file_dir, 'PEDV-rna-cell-2rep.1nt.none.matrix', res)[0] 
ma_sars_rna = extract_matrix(file_dir, 'SARS-cov2-all-rep-rna.1nt.none.matrix', res)[0]
ma_ibv_rna_cell = extract_matrix(file_dir, 'IBV-rna.1nt.none.matrix', res)[0]
ma_pdcov_rna = extract_matrix(file_dir, 'WT_in_Virion.none.1bp.matrix', res)[0]
ma_sars1_rna = extract_matrix(file_dir, 'Sars-cov-24hpi.1nt.none.matrix', res)[0]
ma_mers_rna = extract_matrix(file_dir, 'PRJNA279442-24h_MERS-rnaseq-PRJNA279442-calu3-24h.1nt.none.matrix', res)[0]

cmap = 'afmhot_r' #'PuBu'
# joint_heatmap_hist(ma_pedv_rna_cell, 'PEDV', 10, 5, cmap)
# joint_heatmap_hist(ma_sars_rna, 'SARS-CoV2', 10, 5, cmap)
# joint_heatmap_hist(ma_ibv_rna_cell, 'IBV', 10, 0.25, cmap)
# joint_heatmap_hist(ma_pdcov_rna, 'PDCoV', 10, 5, cmap)

# sns.set_theme(style="ticks")
# pp = PdfPages('G:/OneDrive - webmail.hzau.edu.cn/vRic-seq/Fig-ORF10/geneme-junction-sites-new3.pdf')
# pp.savefig(joint_heatmap_hist_genome(ma_pdcov_rna, dict_gene_sites_pdcov, 'PDCoV', res, 10, cmap))
# pp.savefig(joint_heatmap_hist_genome(ma_ibv_rna_cell, dict_gene_sites_ibv, 'IBV', res, 2.5, cmap))
# pp.savefig(joint_heatmap_hist_genome(ma_pedv_rna_cell, dict_gene_sites_pedv, 'PEDV', res, 17.5, cmap))
# pp.savefig(joint_heatmap_hist_genome(ma_mers_rna, dict_gene_sites_mers, 'MERS-CoV', res, 15, cmap))
# pp.savefig(joint_heatmap_hist_genome(ma_sars_rna, dict_gene_sites_sars, 'SARS-CoV-2', res, 15, cmap))
# pp.savefig(joint_heatmap_hist_genome(ma_sars1_rna, dict_gene_sites_sars1, 'SARS-CoV', res, 8, cmap))
# pp.close()

pedv_utr3 = dict_gene_sites_pedv['UTR_3'][1] - dict_gene_sites_pedv['UTR_3'][0]
sars_utr3 = dict_gene_sites_sars['UTR_3'][1] - dict_gene_sites_sars['UTR_3'][0]
ibv_utr3 = dict_gene_sites_ibv['UTR_3'][1] - dict_gene_sites_ibv['UTR_3'][0]
pdcov_utr3 = dict_gene_sites_pdcov['UTR_3'][1] - dict_gene_sites_pdcov['UTR_3'][0]
sars1_utr3 = dict_gene_sites_sars1['UTR_3'][1] - dict_gene_sites_sars1['UTR_3'][0]
mers_utr3 = dict_gene_sites_mers['UTR_3'][1] - dict_gene_sites_mers['UTR_3'][0]

sns.set_theme(style="whitegrid")
fig = plt.figure(figsize=(8,6))
plt.bar(np.arange(1,7) - 0.155, [pedv_utr3, mers_utr3, sars_utr3, sars1_utr3, ibv_utr3, pdcov_utr3], .31, color = 'C0', alpha = 0.75, label = 'Length(bp)')
plt.xticks([1,2,3,4,5,6],['PEDV', 'MERS-CoV', 'SARS-CoV2', 'SARS-CoV', 'IBV', 'PDCoV'])
ax = plt.gca()
ax.spines['top'].set_visible(False)
ax.set_ylabel('Length(bp)')
ax2 = ax.twinx()
plt.bar(np.arange(1,7) + 0.155, [pedv_utr3/ dict_gene_sites_pedv['UTR_3'][1], mers_utr3/ dict_gene_sites_mers['UTR_3'][1], sars_utr3/ dict_gene_sites_sars['UTR_3'][1]
                    , sars1_utr3/ dict_gene_sites_sars['UTR_3'][1], ibv_utr3/  dict_gene_sites_ibv['UTR_3'][1], pdcov_utr3/  dict_gene_sites_pdcov['UTR_3'][1]], .31
        , color = 'C0', alpha = 0.5, label = 'Proportion')
ax2.spines['top'].set_visible(False)
ax2.set_ylabel('Proportion')
fig.legend(loc='upper center')



