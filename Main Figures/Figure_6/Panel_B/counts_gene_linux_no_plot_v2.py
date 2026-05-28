# -*- coding: utf-8 -*-
"""
Created on Tue Nov  8 11:12:36 2022

@author: DELL
"""
from collections import OrderedDict, Counter
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from sys import argv

pyname, chrs, file_in, file_name = argv 

#dict_gene = OrderedDict([('leader', (0, 100)),

#genes = ['ORF1a', 'ORF1b', 'S', 'ORF3a', 'E', 'M', 'ORF6', 'ORF7a', 'ORF7b', 'ORF8', 'N', 'ORF9b', 'ORF9c', 'ORF10']
#genes = list(dict_gene.keys())[2:-1]

# chrs = 'MK584552'
# file_dir = 'G:/OneDrive - webmail.hzau.edu.cn/vRic-seq/PEDV/nanopore/'
# file_name = 'MK584552.1.raw.splice'
dict_gene = OrderedDict([('leader', (1, 100)),
            ('UTR_5', (101, 265)),
            ('ORF1a', (266,13467)),
            ('ORF1b', (13468,21552)),
            ('S', (21563,25384)),
            ('ORF3a', (25393,26220)),
            ('E', (26245,26472)),
            ('M', (26523,27191)),
            ('ORF6', (27202,27387)),
            ('ORF7a', (27394,27759)),
            ('ORF7b', (27756,27887)),
            ('ORF8',(27894,28259)),
            ('N',(28274,29533)),
            ('ORF10',(29558,29674)),
            ('UTR_3', (29675,29903))])

def is_in(dict_data, in_data):
    out_data = []
    for key in dict_data:
        values = dict_data[key]
        if values[0] <=in_data and values[1] >= in_data:
            out_data.append(key)
    return out_data

def is_between(dict_data, in_data):
    out_data = []
    for i,key in enumerate(dict_data):
        values = dict_data[key]
        
        if i == 0:
            if values[0] <=in_data[0] and values[1] >= in_data[0]:
                #片段起点是否在leather区域内
                out_data.append(key)
        
        else:
            if values[0] >=in_data[0] and values[1] <= in_data[1]:
                #判断片段是否覆盖整个区域
                out_data.append(key)
                
    return out_data

def is_leader(dict_data, in_data):
    out_data = []
    for i,key in enumerate(dict_data):
        values = dict_data[key]
        
        if i == 0:
            if values[0] <=in_data[0] and values[1] >= in_data[0]:
                out_data.append(key)
        
        else:
            if values[0] >=in_data[0] and values[1] <= in_data[1]:
                out_data.append(key)
                
    return out_data

def junction(dict_data, in_data):
    out_data = []
    for key in dict_data:
        #判断片段起点在那个区域内
        values = dict_data[key]
        if values[0] <=in_data[0] and values[1] >= in_data[0]:
            out_data.append(key)
        
    if len(out_data) == 0:
        #片段起点没在区域内就是gap区域    
        out_data.append('gap')
        
    for key in dict_data:
        #判断片段终点在那个区域内
        values = dict_data[key]        
        if values[0] <=in_data[1] and values[1] >= in_data[1]:
            out_data.append(key)
            
    if len(out_data) == 1:    
            out_data.append('gap')
            
    return '-'.join(out_data)

def delete_list(list_in, list_delet):
    for x in list_delet:
        if x in list_in:
            list_in.remove(x)
    return list_in
    

files_in = open(file_in)
files_out2 = open(file_name + '-2.bed', 'w')
files_out3 = open(file_name + '-3.bed', 'w')
max_junc = 0
for line in files_in:
    lines = line.split('\t') 
    reads_name = lines.pop(0)
    lens = len(lines)
    iters = int(lens/2) #片段对数
    names = [] #起点终点
    trscpt = [] #全覆盖区域
    trscpt_junc = [] #起点终点所在区域
    if lens % 2 != 0:
        print('die out')
        break
    for i in range(iters):
        start = int(lines[i*2])
        end = int(lines[i*2+1])
        if i == (iters-1) :
            end = dict_gene['UTR_3'][1] #把最后一个片段的终点变为UTR3的终点
        gene_name = is_between(dict_gene, [start, end])
        #判断片段是否在leather中以及判断片段是否覆盖整个区域
        trscpt.extend(gene_name)
        trscpt_junc.append(junction(dict_gene,[start, end]))
        #判断片段起点在那个区域内，判断片段终点在那个区域内
        names.extend( [str(start), str(end)]) 
        
    trscpt = delete_list(trscpt, ['leader', 'UTR_5'])
        
    if len(trscpt) == 0 :
        continue
    juncs = iters - 1 #跳了多少次
    if juncs > max_junc :
        max_junc = juncs

    values = '->'.join(trscpt_junc)
                
    # if true_gen == 'ORF9b':
    if len(names) >= 5:
        #跳2次以上
        files_out3.writelines(reads_name+'\t'+'\t'.join(names)+'\t'+values+'\n')
    else :
        #跳2次
        files_out2.writelines(reads_name+'\t'+'\t'.join(names)+'\t'+values+'\n')
    # if len(bed) >= 7:
    #     files_out.writelines('\t'.join(bed)+'\n')
files_out2.close()
files_out3.close()
files_in.close()
#print(Counter(transcpt_out[('ORF10', 1)]))