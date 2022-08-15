import os
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as st
import statsmodels
#from statsmodels.stats.multicomp import pairwise_tukeyhsd
#from statsmodels.stats.libqsturng import psturng
from pandas import DataFrame, Categorical, Series
from statsmodels.sandbox.stats.multicomp import multipletests
import seaborn as sns
import itertools as it
from itertools import compress
from statsmodels.compat.python import range
from statsmodels.compat.collections import OrderedDict
#import time
import math
from scipy.stats import rankdata
import pylab
import statsmodels.graphics.gofplots as sm
import pickle
import pandas as pd
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import operator
from scipy.interpolate import make_interp_spline
from scipy.stats import fisher_exact


def take(n, iterable):
    "Return first n items of the iterable as a list"
    return list(islice(iterable, n))

def posthoc_conover(a, val_col = None, group_col = None, p_adjust = None, sort = True):

    '''Post-hoc pairwise test for multiple comparisons of mean rank sums
    (Conover's test). May be used after Kruskal-Wallis one-way analysis of
    variance by ranks to do pairwise comparisons.
        Parameters
        ----------
        a : array_like or pandas DataFrame object
            An array, any object exposing the array interface or a pandas DataFrame.
            Array must be two-dimensional. Second dimension may vary,
            i.e. groups may have different lengths.
        val_col : str, optional
            Must be specified if `a` is a pandas DataFrame object.
            Name of the column that contains values.
        group_col : str, optional
            Must be specified if `a` is a pandas DataFrame object.
            Name of the column that contains group names.
        p_adjust : str, optional
            Method for adjusting p values. See statsmodels.sandbox.stats.multicomp
            for details. Available methods are:
                'bonferroni' : one-step correction
                'sidak' : one-step correction
                'holm-sidak' : step-down method using Sidak adjustments
                'holm' : step-down method using Bonferroni adjustments
                'simes-hochberg' : step-up method  (independent)
                'hommel' : closed method based on Simes tests (non-negative)
                'fdr_bh' : Benjamini/Hochberg  (non-negative)
                'fdr_by' : Benjamini/Yekutieli (negative)
                'fdr_tsbh' : two stage fdr correction (non-negative)
                'fdr_tsbky' : two stage fdr correction (non-negative)
        sort : bool, optional
            Specifies whether to sort DataFrame by group_col or not. Recommended
            unless you sort your data manually.
        Returns
        -------
        Numpy ndarray if `a` is an array-like object else pandas DataFrame of p values.
        Notes
        -----
        A tie correction are employed according to Conover (1979).
        References
        ----------
        W. J. Conover and R. L. Iman (1979), On multiple-comparisons procedures, Tech.
        Rep. LA-7677-MS, Los Alamos Scientific Laboratory.
        Examples
        --------
    '''

    def compare_conover(i, j):
        diff = np.abs(x_ranks_avg[i] - x_ranks_avg[j])
        B = (1. / x_lens[i] + 1. / x_lens[j])
        D = (x_len_overall - 1. - H) / (x_len_overall - x_len)
        t_value = diff / np.sqrt(S2 * B * D)
        p_value = 2. * st.t.sf(np.abs(t_value), df = x_len_overall - x_len)
        return p_value

    def get_ties(x):
        x_sorted = np.array(np.sort(x))
        tie_sum = 0
        pos = 0
        while pos < x_len_overall:
            n_ties = len(x_sorted[x_sorted == x_sorted[pos]])
            pos = pos + n_ties
            if n_ties > 1:
                tie_sum += n_ties ** 3. - n_ties
        c = np.min([1., 1. - tie_sum / (x_len_overall ** 3. - x_len_overall)])
        return c

    if isinstance(a, DataFrame):
        x = a.copy()
        if not sort:
            x[group_col] = Categorical(x[group_col], categories=x[group_col].unique(), ordered=True)
        x.sort_values(by=[group_col, val_col], ascending=True, inplace=True)
        x_groups_unique = np.asarray(x[group_col].unique())
        x_len = x_groups_unique.size
        x_lens = x.groupby(by=group_col)[val_col].count().values
        x_flat = x[val_col].values
        x_lens_cumsum = np.insert(np.cumsum(x_lens), 0, 0)[:-1]
        x_grouped = np.array([x_flat[j:j + x_lens[i]] for i, j in enumerate(x_lens_cumsum)])

    else:
        x = np.array(a)
        x_grouped = np.array([np.asarray(a)[~np.isnan(a)] for a in x])
        x_flat = np.concatenate(x_grouped)
        x_len = len(x_grouped)
        x_lens = np.asarray([len(a) for a in x_grouped])
        x_lens_cumsum = np.insert(np.cumsum(x_lens), 0, 0)[:-1]

    x_len_overall = len(x_flat)

    if any(x_lens == 0):
        raise ValueError("All groups must contain data")

    x_ranks = st.rankdata(x_flat)
    x_ranks_grouped = np.array([x_ranks[j:j + x_lens[i]] for i, j in enumerate(x_lens_cumsum)])
    x_ranks_avg = [np.mean(z) for z in x_ranks_grouped]
    x_ties = get_ties(x_ranks) #st.tiecorrect(x_ranks)

    H = st.kruskal(*x_grouped)[0]

    if x_ties == 1:
        S2 = x_len_overall * (x_len_overall + 1.) / 12.
    else:
        S2 = (1. / (x_len_overall - 1.)) * (np.sum(x_ranks ** 2.) - (x_len_overall * (((x_len_overall + 1.)**2.) / 4.)))

    vs = np.zeros((x_len, x_len), dtype=np.float)
    tri_upper = np.triu_indices(vs.shape[0], 1)
    tri_lower = np.tril_indices(vs.shape[0], -1)
    vs[:,:] = 0

    combs = it.combinations(range(x_len), 2)

    for i,j in combs:
        vs[i, j] = compare_conover(i, j)

    if p_adjust:
        vs[tri_upper] = multipletests(vs[tri_upper], method = p_adjust)[1]
    vs[tri_lower] = vs.T[tri_lower]
    np.fill_diagonal(vs, -1)

    if isinstance(x, DataFrame):
        return DataFrame(vs, index=x_groups_unique, columns=x_groups_unique)
    else:
        return vs




a1 = 2

directory = sys.argv[1]

"""
#This block makes a dictionary with each gene having binding positions labeled {gene:[1,0,1,0,0...],gene:...}
"""

domain_positions = pickle.load(open(directory+"full_domain_positions.pickle","rb"))


f = open('/Genomics/grid/users/sbking/run_pipeline/primate_bind_ens_human_proteins_1.0_1.0.fasta' , 'r')
bin_dic = {}
counter = 0

for line in f:
	if '>' in line:
		gene = line.rstrip()
		bin_dic[gene] = []
		dom_pos_count = 1 #it is important that you get the indexing correct here, hmmer output is 1 indexed
	else:
		linesplit = line.rstrip().split(',')
		for each in linesplit:
			if str(each[0]) == 'X':
				dom_pos_count += 1
				continue
			try:
				if dom_pos_count in domain_positions[gene]:
					if int(each[-1]) == 0:
						bin_dic[gene].append(a1) #switch this between 0 and 2 to change taking into account domains or not (really no reason to turn this off)
					else:
						bin_dic[gene].append(int(each[-1]))
				else:
					bin_dic[gene].append(int(each[-1]))
				dom_pos_count += 1
			except:
				bin_dic[gene].append(int(each[-1]))
				dom_pos_count +=1
f.close()

os.chdir(directory+'final_orthogroups')

quick_fils = os.listdir('.')

ortho_gene_dic = {}
reverse_ortho_gene_dic = {}

for fil in quick_fils:
	os.chdir(directory+'final_orthogroups/'+fil)
	new_quick_fils = os.listdir('.')
	for each in new_quick_fils:
		if 'cat_' in each:
			fg = open(each,'r')
			for line in fg:
				if '>9606:' in line:
					q_gn = line.rstrip().split(' ')[0]
					ortho_gene_dic[fil] = q_gn
					reverse_ortho_gene_dic[q_gn.split('_')[1]] = fil

			fg.close()


os.chdir(directory+'final_orthogroups')

fils = os.listdir('.')
pvals = []
fil_val = []
tot_rank_dic = {}
per_prot_vals = []
tot_score_list = []
	
fil_count = 0
r_count = 0

fab = open(directory+'background_gene_list.txt','w')


for fil in fils:
	gn = ortho_gene_dic[fil]
	fab.write((gn.split('_')[1])+'\n')
	per_gene_rank_dic = {}
	newfile = fil.rstrip()
	bind = []
	nonbind = []
	os.chdir('./'+fil)
	newfils = os.listdir('.')
	fil_count += 1
	for each in newfils:
		if each == ('beat10.mafft'):
			total_res_count = []
			res_count = 0
			z = open(each, 'r')
			for li in z:
				res_count = 0
				if ('>' not in li) and (li != '\n') and (li != '') and (li != ' '):
					li = li.rstrip()
					for residue in li:
						res_count = res_count + 1
					total_res_count.append(res_count)
			z.close()	
	g = open('r4sOrig.res', 'r')
	per_prot_vals = []
	for line in g:
		if '[' in line:
			linesplit = line.rstrip().split(' ')
			linesplit = [x for x in linesplit if x!= '']
			rs = int(linesplit[0]) - 1        #this line makes it 0 indexed
			tot_score_list.append(float(linesplit[2]))
			per_prot_vals.append(float(linesplit[2]))
			per_gene_rank_dic[gn+'_'+repr(rs)] = [float(linesplit[2]), bin_dic[gn][rs]]
			r_count += 1
	
	
	g.close()



	se = open('SE.txt_last', 'r')

	for line in se:
		if "#" in line:
			continue
		linesplit = line.rstrip().split('\t')
		rs = int(linesplit[0]) - 1
		cons = float(linesplit[2])
		per_gene_rank_dic[gn+'_'+repr(rs)].append(cons)


	se.close()


	tot_rank_dic.update(per_gene_rank_dic) #r4s score, binding (binary), se score


	os.chdir('../')


dump_dic = open(directory+'domain_orig_r4s_dic.pickle', 'wb')

pickle.dump(tot_rank_dic,dump_dic)

fab.close()

binding_conserved = 0
binding_total = 0
nonbinding_conserved = 0
nonbinding_total = 0
domain_non_conserved = 0
domain_non_total = 0


variantz = pickle.load(open(directory+"cutoff_variants_dic_of_chroms.p", "rb"))

enrichment_dic = {}

for gene in tot_rank_dic:
	enrichment_dic[gene.split('_')[1]] = [0,0,0,0,0]   #total variable, total, binding variable, binding total, fold enrichment


for gene in tot_rank_dic:
	if gene.split('_')[1] in reverse_ortho_gene_dic:
		pass
	else:
		continue
	if tot_rank_dic[gene][1] == 0:
		if tot_rank_dic[gene][2] < 0.005:
			nonbinding_conserved += 1
		nonbinding_total += 1
	if tot_rank_dic[gene][1] == 1:
		if tot_rank_dic[gene][2] < 0.005:
			binding_conserved += 1
		binding_total += 1
	if tot_rank_dic[gene][1] == 2:
		if tot_rank_dic[gene][2] < 0.005:
			domain_non_conserved += 1
		else:
			if gene.split('_')[1] not in variantz:
				enrichment_dic[gene.split('_')[1]][0] += 1
			else:
				if gene.split('_')[-1] not in variantz[gene.split('_')[1]]:
					enrichment_dic[gene.split('_')[1]][0] += 1
		enrichment_dic[gene.split('_')[1]][1] += 1
		domain_non_total += 1



with open(directory+'total_dict.pickle', 'wb') as handle:
	pickle.dump(tot_rank_dic, handle)



domain_binding_group = []
domain_nonbinding_group = []
nondomain_nonbinding_group = []

domain_binding_group_SE = []
domain_nonbinding_group_SE = []
nondomain_nonbinding_group_SE = []


for each_thing in tot_rank_dic:
	if tot_rank_dic[each_thing][1]==1:
		domain_binding_group.append(tot_rank_dic[each_thing][0])
		domain_binding_group_SE.append(tot_rank_dic[each_thing][2])
	if tot_rank_dic[each_thing][1]==0:
		nondomain_nonbinding_group.append(tot_rank_dic[each_thing][0])
		nondomain_nonbinding_group_SE.append(tot_rank_dic[each_thing][2])
	if tot_rank_dic[each_thing][1]==2:
		domain_nonbinding_group.append(tot_rank_dic[each_thing][0])
		domain_nonbinding_group_SE.append(tot_rank_dic[each_thing][2])


type_list = ['SM','DRUGLIKE','METABOLITE','ION','PEPTIDE','NUCACID','DNA','DNABASE','DNABACKBONE','RNA','RNABASE','RNABACKBONE']


dic_of_types = {}

global_number_binding = 0
global_number_conserved = 0

listofusedsites = []

paper_print_dic = {}


num_ortho_set = set()

for thing in type_list:
	annoying = True
	linecount = 0
	counts_dics = {}
	bin_poses_dic = {}
	dic_of_types[thing] = [[],[],[]] #r4s, binding, SE
	number_conserved = 0
	number_binding = 0
	ff = open(directory+'without_scores/primate_'+thing+'_list_1.0_1.0', 'r')
	zzz = open(directory+'rapid_types/'+thing+'.txt', 'w')
	for line in ff:
		line = line.rstrip()
		if line not in tot_rank_dic:
			continue
		if tot_rank_dic[line][1] == 1 or annoying:
			enrichment_dic[line.split('_')[1]][3] += 1
			enrichment_dic[line.split('_')[1]][1] += 1
			linecount += 1
			number_binding += 1
			if tot_rank_dic[line][0] > 2.5: #this number here is the evolutionary rate, next part does rapid sites
				if line.split('_')[1] not in bin_poses_dic:
					bin_poses_dic[line.split('_')[1]] = [str(int(line.split('_')[-1])+1)]
				else:
					bin_poses_dic[line.split('_')[1]].append(str(int(line.split('_')[-1])+1))
				if line.split('_')[1] not in counts_dics:
					counts_dics[line.split('_')[1]] = 1
				else:
					counts_dics[line.split('_')[1]] = counts_dics[line.split('_')[1]] + 1
			if tot_rank_dic[line][2] > 0.005: #this is to only grab sites that are varying, you can turn this off if you have the need to
				num_ortho_set.add(line.split('_')[1])
				if line.split('_')[3] not in variantz:
					enrichment_dic[line.split('_')[1]][2] += 1
					enrichment_dic[line.split('_')[1]][0] += 1
				else:
					if line.split('_')[3] not in variantz[gene.split('_')[1]]:
						enrichment_dic[line.split('_')[1]][2] += 1
						enrichment_dic[line.split('_')[1]][0] += 1
			dic_of_types[thing][0].append((tot_rank_dic[line][0]))
			dic_of_types[thing][1].append((tot_rank_dic[line][1]))
			dic_of_types[thing][2].append((tot_rank_dic[line][2]))
			if tot_rank_dic[line][2] < 0.005:
				number_conserved += 1
			if line in listofusedsites:
				pass
			else:
				if tot_rank_dic[line][2] < 0.005:
					global_number_conserved += 1
				global_number_binding += 1
				listofusedsites.append(line)
			
	paper_print_dic[thing] = (linecount, (number_binding), float(number_conserved)/number_binding)


	ff.close()


	for eachh in counts_dics:
		out_str = " "
		zzz.write(eachh + '\t' + repr(counts_dics[eachh]) + '\t' + reverse_ortho_gene_dic[eachh] + '\t'+out_str.join(bin_poses_dic[eachh])+'\n') 
			
	zzz.close()




f_gene = open(directory+"genes_for_GO.txt", 'w')

for each_gene in num_ortho_set:
	f_gene.write(each_gene + '\n')


f_gene.close()



print ('Number of ligand-binding sites: ', global_number_binding)
print ('Number of total other sites (both outside and in domains): ', (nonbinding_total+domain_non_total))
print ('Number of other sites in domains: ', domain_non_total)

print ('Number of primate ortho-groups: ',fil_count)


print ('Percent of ligand-binding sites conserved: ', 100*float(global_number_conserved/global_number_binding))
print ('Percent of other sites outside of domains conserved: ',100*(float((domain_non_conserved+nonbinding_conserved))/(domain_non_total+nonbinding_total)))
print ('Percent of other sites in domains conserved: ', 100*(float(domain_non_conserved)/domain_non_total))

print ('Ligand-binding vs other sites (both outside and in domains) percent conservation Fishers exact test: ', fisher_exact(np.array([[global_number_conserved, global_number_binding],[(domain_non_conserved+nonbinding_conserved),(domain_non_total+nonbinding_total)]])))

print ('Other domain sites vs other sites (both outside and in domains) percent conservation Fishers exact test: ', fisher_exact(np.array([[domain_non_conserved, domain_non_total],[(domain_non_conserved+nonbinding_conserved),(domain_non_total+nonbinding_total)]])))

print ('Ligand-binding vs other domain sites percent conservation Fishers exact test: ', fisher_exact(np.array([[global_number_conserved, global_number_binding],[(domain_non_conserved),(domain_non_total)]])))


names_dic = {'SM':"Small Molecule", 'ION':"Ion", 'DNA':"DNA", 'RNA':"RNA", 'PEPTIDE':"Peptide"}

for each in paper_print_dic:

	if each in names_dic:
		print ('% variable sites: ',names_dic[each], 100-(100*float(paper_print_dic[each][2])))

print('Ligand-binding vs other sites evolutionary rate means and Mann Whitney-U test: ',np.mean(domain_binding_group), np.mean(domain_nonbinding_group), st.mannwhitneyu(domain_binding_group,domain_nonbinding_group))

print('Ligand-binding vs other sites shannon entropy means and Mann Whitney-U test: ',np.mean(domain_binding_group_SE), np.mean(domain_nonbinding_group_SE), st.mannwhitneyu(domain_binding_group_SE,domain_nonbinding_group_SE))



names_dic = {'SM':"Small Molecule", 'ION':"Ion", 'DNA':"DNA", 'RNA':"RNA", 'PEPTIDE':"Peptide"}

for each in paper_print_dic:
	
	if each in names_dic:
		print ('Number of sites and % that varies :',names_dic[each], paper_print_dic[each][1], 100*float(paper_print_dic[each][2]))



enr_count = 0

enr_list = []

f = open(directory+"enriched_for_GO.txt", 'w')


for each in enrichment_dic:
	total_var = enrichment_dic[each][0]
	total = enrichment_dic[each][1] 
	bin_var = enrichment_dic[each][2]
	bin_total = enrichment_dic[each][3] 

	if total == 0 or bin_total == 0:
		continue
	
	tot_enrich = float(total_var)/total
	bin_enrich = float(bin_var)/bin_total

	if tot_enrich == 0:
		if bin_enrich != 0:
			pass
		continue
	if float(bin_var) < 2:
		continue

	enrichment_dic[each][4] = float(bin_enrich)/tot_enrich
	enr_list.append(float(bin_enrich)/tot_enrich)

	if enrichment_dic[each][4] >= 2:
		enr_count += 1
		f.write(each+'\n')


f.close()


print ('Number of ortho-groups with 2-fold+ enriched variable sites: ', enr_count)

with open(directory+'dic_of_types.pickle', 'wb') as handle:
	pickle.dump(dic_of_types, handle, protocol=2)





krus_res = st.kruskal(dic_of_types['SM'][0],dic_of_types['DRUGLIKE'][0],dic_of_types['METABOLITE'][0],dic_of_types['ION'][0],dic_of_types['PEPTIDE'][0],dic_of_types['NUCACID'][0],dic_of_types['DNA'][0],dic_of_types['DNABASE'][0],dic_of_types['DNABACKBONE'][0],dic_of_types['RNA'][0],dic_of_types['RNABASE'][0],dic_of_types['RNABACKBONE'][0])



krus_p = krus_res[1]

tot_arr_all_bind = [dic_of_types['SM'][0],dic_of_types['DRUGLIKE'][0],dic_of_types['METABOLITE'][0],dic_of_types['ION'][0],dic_of_types['PEPTIDE'][0],dic_of_types['NUCACID'][0],dic_of_types['DNA'][0],dic_of_types['DNABASE'][0],dic_of_types['DNABACKBONE'][0],dic_of_types['RNA'][0],dic_of_types['RNABASE'][0],dic_of_types['RNABACKBONE'][0]]

if krus_p < 0.001:
	conover_res = posthoc_conover(tot_arr_all_bind, p_adjust = 'bonferroni')



means = [np.mean(dic_of_types['SM'][0]),np.mean(dic_of_types['DRUGLIKE'][0]),np.mean(dic_of_types['METABOLITE'][0]),np.mean(dic_of_types['ION'][0]),np.mean(dic_of_types['PEPTIDE'][0]),np.mean(dic_of_types['NUCACID'][0]),np.mean(dic_of_types['DNA'][0]),np.mean(dic_of_types['DNABASE'][0]),np.mean(dic_of_types['DNABACKBONE'][0]),np.mean(dic_of_types['RNA'][0]),np.mean(dic_of_types['RNABASE'][0]),np.mean(dic_of_types['RNABACKBONE'][0])]


means_se = [np.mean(dic_of_types['SM'][2]),np.mean(dic_of_types['DRUGLIKE'][2]),np.mean(dic_of_types['METABOLITE'][2]),np.mean(dic_of_types['ION'][2]),np.mean(dic_of_types['PEPTIDE'][2]),np.mean(dic_of_types['NUCACID'][2]),np.mean(dic_of_types['DNA'][2]),np.mean(dic_of_types['DNABASE'][2]),np.mean(dic_of_types['DNABACKBONE'][2]),np.mean(dic_of_types['RNA'][2]),np.mean(dic_of_types['RNABASE'][2]),np.mean(dic_of_types['RNABACKBONE'][2])]


medians = [np.median(dic_of_types['SM'][0]),np.median(dic_of_types['DRUGLIKE'][0]),np.median(dic_of_types['METABOLITE'][0]),np.median(dic_of_types['ION'][0]),np.median(dic_of_types['PEPTIDE'][0]),np.median(dic_of_types['NUCACID'][0]),np.median(dic_of_types['DNA'][0]),np.median(dic_of_types['DNABASE'][0]),np.median(dic_of_types['DNABACKBONE'][0]),np.median(dic_of_types['RNA'][0]),np.median(dic_of_types['RNABASE'][0]),np.median(dic_of_types['RNABACKBONE'][0])]


medians_se = [np.median(dic_of_types['SM'][2]),np.median(dic_of_types['DRUGLIKE'][2]),np.median(dic_of_types['METABOLITE'][2]),np.median(dic_of_types['ION'][2]),np.median(dic_of_types['PEPTIDE'][2]),np.median(dic_of_types['NUCACID'][2]),np.median(dic_of_types['DNA'][2]),np.median(dic_of_types['DNABASE'][2]),np.median(dic_of_types['DNABACKBONE'][2]),np.median(dic_of_types['RNA'][2]),np.median(dic_of_types['RNABASE'][2]),np.median(dic_of_types['RNABACKBONE'][2])]


counti = 0
for i in conover_res:
	countj = 0
	for j in i:
		if j >= 0.001:
			conover_res[counti][countj] = 0
			countj = countj + 1
			continue
		if j < 0.001:
			if means[counti] >= means[countj]:  
				conover_res[counti][countj] = 0
				countj = countj + 1
				continue
		countj = countj + 1
	counti = counti + 1


counti = 0
for i in conover_res:
        countj = 0
        for j in i:
                if float(j) > 0:
                        conover_res[counti][countj] = -1*(np.log10(conover_res[counti][countj]))
                countj += 1
        counti += 1




