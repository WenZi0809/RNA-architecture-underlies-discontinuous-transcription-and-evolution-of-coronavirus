# -*- coding: utf-8 -*-
"""
Created on Tue Jun 20 10:48:27 2023

@author: Administrator
"""

from collections import OrderedDict
import matplotlib.pyplot as plt
import numpy as np

file_dir = 'G:/OneDrive - webmail.hzau.edu.cn/vRic-seq/data/genome/'
color = ['lightgreen','lightgreen', 'pink', 'lightblue', 'pink', 'lightblue', 'pink', 'lightblue', 'pink', 'lightblue','pink', 'lightblue','pink', 'lightblue']

gene_sites_sars = np.loadtxt(file_dir + 'NC_045512.2.bed', usecols=[3,1,2], 
                        dtype = [('gene', '<U6'),('S', '<i4'),('E', '<i4')], skiprows=1)
gene_sites_pedv = np.loadtxt(file_dir + 'MK584552.1.bed', usecols=[3,1,2], 
                        dtype = [('gene', '<U6'),('S', '<i4'),('E', '<i4')], skiprows=1)
gene_sites_pdcov_kv = np.loadtxt(file_dir + 'KY363867.1.bed', usecols=[3,1,2], 
                        dtype = [('gene', '<U6'),('S', '<i4'),('E', '<i4')], skiprows=1)
gene_sites_pdcov_kt = np.loadtxt(file_dir + 'KT336560.1.bed', usecols=[3,1,2], 
                        dtype = [('gene', '<U6'),('S', '<i4'),('E', '<i4')], skiprows=1)
gene_sites_ibv = np.loadtxt(file_dir + 'MK878536.1.bed', usecols=[3,1,2], 
                        dtype = [('gene', '<U6'),('S', '<i4'),('E', '<i4')], skiprows=1)
gene_sites_sars1 = np.loadtxt(file_dir + 'NC_004718.3.bed', usecols=[3,1,2], 
                        dtype = [('gene', '<U6'),('S', '<i4'),('E', '<i4')], skiprows=1)
gene_sites_mers = np.loadtxt(file_dir + 'NC_019843.3.bed', usecols=[3,1,2], 
                        dtype = [('gene', '<U6'),('S', '<i4'),('E', '<i4')], skiprows=1)

dict_gene_sites_sars = OrderedDict()
for gene in gene_sites_sars:
    dict_gene_sites_sars[gene[0]] = (gene[1], gene[2]) 

dict_gene_sites_pedv = OrderedDict()
for gene in gene_sites_pedv:
    dict_gene_sites_pedv[gene[0]] = (gene[1], gene[2]) 
    
dict_gene_sites_pdcov = OrderedDict()
for gene in gene_sites_pdcov_kt:
    dict_gene_sites_pdcov[gene[0]] = (gene[1], gene[2]) 

dict_gene_sites_ibv = OrderedDict()
for gene in gene_sites_ibv:
    dict_gene_sites_ibv[gene[0]] = (gene[1], gene[2]) 

dict_gene_sites_sars1 = OrderedDict()
for gene in gene_sites_sars1:
    dict_gene_sites_sars1[gene[0]] = (gene[1], gene[2])  

dict_gene_sites_mers = OrderedDict()
for gene in gene_sites_mers:
    dict_gene_sites_mers[gene[0]] = (gene[1], gene[2])  

