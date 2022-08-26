#######################
# GENERAL INFORMATION #
#######################
1. Title of Dataset: Code and data from "Primate protein-ligand interfaces exhibit significant conservation and unveil human-specific evolutionary drivers"

2. Author Information
	A. Principal Investigator Contact Information
		Name: Mona Singh
		Institution: Princeton University, Department of Computer Science
		Address: Princeton, NJ 08540, USA
		Email: mona@cs.princeton.edu

	B. Associate or Co-investigator Contact Information
		Name: Sean King
		Institution: Princeton University, Department of Molecular Biology
		Address: Princeton NJ 08540, USA
		Email: sbking@princeton.edu

3. Date of data collection and code production (2018-03-01 - 2022-06-01)

4. Geographic location of code production: Princeton, NJ

5. Information about funding sources that supported the collection of the data: NIH (grant R01-GM76275) to Mona Singh and NIH (grant 5T32GM007388) to Princeton University Department of Molecular Biology

##############################
# SHARING/ACCESS INFORMATION #
##############################

1. Licenses/restrictions placed on the data: MIT License

Copyright (c) [2022] [Sean Burak King]

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

2. Links to publications that cite or use the data: King S., and Singh M. Primate protein-ligand interfaces exhibit significant conservation and unveil human-specific evolutionary drivers. 2022. In preparation.

3. Links to other publicly accessible locations of the data: The three accessible locations for the data include: 1. The manuscript, 2. The supplementary information associated with the manuscript, and 3. This repository

4. Updates: This final version of this database now includes scripts and data that calculate all results found within the manuscript

5. Links/relationships to ancillary data sets: All required and available data is included in this dataspace repository

6. Was data derived from another source? yes/no: No

7. Recommended citation for this dataset:
	Dataset: King S., and Singh M. Code and data from "Primate protein-ligand interfaces exhibit significant conservation and unveil human-specific evolutionary drivers". https://github.com/Singh-Lab/The-evolutionary-landscape-of-primate-protein-ligand-interfaces. Deposited 2022.
	Publication: King S., and Singh M. Primate protein-ligand interfaces exhibit significant conservation and unveil human-specific evolutionary drivers. 2022. In preparation.


########################
# DATA & FILE OVERVIEW #
########################

File List:

	
		
		README.txt - this document
		quantify_conservation.py - the instructions on how to run this script are at the bottom of this document
		final_orthogroups, (need to download from included dropbox link)
			-This file contains 10 sub-directories that you must extract into this directory before running the script
			-This will produce 4197 total subfiles for each orthogroup used in the manuscript
			-Each file contains the final alignment used in quantify_conservation.py
			-The raw rate4site (r4sOrig.res) and shannon entropy (SE.txt_last) output is also contained in these files       
        
		background_gene_list.txt - Custom background of genes used in this analysis for unbiased GO search, written by quantify_conservation.py
		primate_bind_ens_human_proteins_1.0_1.0.fasta - per amino-acid binding data for domain sites of each protein used in the analysis
		cutoff_variants_dic_of_chroms.p - dictionary of all human variants under the rare variant cutoff as described in the manuscript
		domain_orig_r4s_dic.pickle.tar.gz - per amino-acid binding data for all domain sites of each protein used in the analysis, written by quantify_conservation.py (you need to unzip this one if you want to use it instead of waiting for the script to generate it)
		total_dict.pickle.tar.gz - per amino-acid binding data for all sites of each protein used in the analysis, written by quantify_conservation.py (you need to unzip this one if you want to use it instead of waiting for the script to generate it)
		without_scores - this file contains per-site list of all ligand binding sites on a per-ligand basis
		rapid_types - list of rapidly evolving sites as described by the manuscript, per orthogroup, segemented by ligand type
		genes_for_GO.txt - list of genes for GO analysis as described in the manuscript, written by quantify_conservation.py
		enriched_for_GO.txt - list of genes with 2x variablity enriched binding sites as described in the manuscript, written by quantify_conservation.py
		dic_of_types.pickle - dictionary of binding sites segmented by the type of ligand they bind, writted by quantify_conservation.py


        
########################
#     INSTRUCTIONS     #
########################


###############################################################################
#                                                                             #
# download the repository to your local file system                           #
# unzip all tar.xz files in the directory, keeping the same name              #
# this master script is provided as a python script for automated running     #
# the only variable you need is the path to your directory                    #
#                                                                             #
# enter the repository and type the command in the line below                 #
# run the command: python3.6.4 quantify_conservation.py /path/to/directory/   #
# to run this script as-is you only need python version 3.6.4                 #                                                                               
#                                                                             #
# this will provide all of the descriptive output in the main manuscript      #
# running this script as-is will give you the exact results in the manuscript #
# the script takes 10-20 minutes to run depending on the platform             #
#                                                                             #
###############################################################################