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

3. Links to other publicly accessible locations of the data: The three accessible locations for the data include: 1. The manuscript, 2. The supplementary information associated with the manuscript, and 3. This dataspace FAIR repository

4. Updates: This final version of this database now includes scripts and data that calculate all results found within the manuscript

5. Links/relationships to ancillary data sets: All required and available data is included in this dataspace repository

6. Was data derived from another source? yes/no: No

7. Recommended citation for this dataset:
	Dataset: King S., and Singh M. Code and data from "Primate protein-ligand interfaces exhibit significant conservation and unveil human-specific evolutionary drivers". Princeton DataSpace. Deposited 2022.
	Publication: King S., and Singh M. Primate protein-ligand interfaces exhibit significant conservation and unveil human-specific evolutionary drivers. 2022. In preparation.

PLEASE UPDATE EVERYTHING BELOW THIS LINE______________________________________________________________________

########################
# DATA & FILE OVERVIEW #
########################

1. File List:

	
		README.txt - this document
		master_script.sh - the instructions on how to run this script to calculate the results in the manuscript are located in the box at the top within the script
		
    PDB Structures
    
		cocomplex_A_6u7e.pdb.gz - PDB co-complex of APN and spike
		cocomplex_S_6m0j.pdb.gz - PDB co-complex of ACE2 and spike       
        
    High Resolution Manuscript Figures
    
		Figure 1.tif - High quality version of figure 1 from the manuscript
		Figure 2.tif - High quality version of figure 2 from the manuscript
        
    Supplementary Information
    
		Supplementary Figure 1.pdf - Example of codon with selection
		Supplementary Figure 2.pdf - NL63 ACE2 sites that bind spike protein
		Supplementary Table 1.xlsx - Total list of protein sequences
		Supplementary Table 2.xlsx - Variation and selection in spike protein
		Supplementary Table 3.xlsx - Reproduction of table 1 in the manuscript without controlling for surface extracellular sites
		Supplementary Table 4.xlsx - Selection within each species group without controlling for surface extracellular sites
		Supplementary Table 5.xlsx - Variation and selection within ACE2
		Supplementary Table 6.xlsx - Variation and selection in APN
		Supplementary Table 7.xlsx - Query and taxid exclusion list

		
	scripts

        apn_r4s_calc.py - calculates z score normalized evolutionary rate estimates for protein group in filename
        apn_se_calc.py - calculates sitewise shannon entropy for protein group in filename
        datamonkey_analysis.py - calculates fold enrichments of pervasive and episodic selection for all groups and a hypergeometric test
        other_minus_r4s_calc.py- calculates z score normalized evolutionary rate estimates for protein group in filename
        other_minus_se_calc.py - calculates sitewise shannon entropy for protein group in filename
        other_r4s_calc.py- calculates z score normalized evolutionary rate estimates for protein group in filename
        other_se_calc.py - calculates sitewise shannon entropy for protein group in filename
        primate_r4s_calc.py- calculates z score normalized evolutionary rate estimates for protein group in filename
        primate_se_calc.py - calculates sitewise shannon entropy for protein group in filename
        rodent_r4s_calc.py- calculates z score normalized evolutionary rate estimates for protein group in filename
        rodent_se_calc.py - calculates sitewise shannon entropy for protein group in filename
        spike_r4s_calc.py- calculates z score normalized evolutionary rate estimates for protein group in filename
        spike_se_calc.py - calculates sitewise shannon entropy for protein group in filename
        total_minus_r4s_calc.py- calculates z score normalized evolutionary rate estimates for protein group in filename
        total_minus_se_calc.py - calculates sitewise shannon entropy for protein group in filename
        total_r4s_calc.py- calculates z score normalized evolutionary rate estimates for protein group in filename
        total_se_calc.py - calculates sitewise shannon entropy for protein group in filename


	data

		total_minus_SE_out - all species minus known SARS hosts shannon entropy output
		total_minus_r4sOrig.res - all species minus known SARS hosts rate4site raw output (this script z score normalizes to match .res output)
		total_SE_out - all species shannon entropy output
		total_r4sOrig.res - all species rate4site raw output (this script z score normalizes to match .res output)
		other_minus_SE_out - "other mammals" species group minus known SARS hosts shannon entropy output
		other_minus_r4sOrig.res - "other mammals" species group minus known SARS hosts rate4site raw output (this script z score normalizes to match .res output)
		other_SE_out - "other mammals" species group shannon entropy output
		other_r4sOrig.res - "other mammals" species group rate4site raw output (this script z score normalizes to match .res output)
		primate_SE_out - primate shannon entropy output
		primate_r4sOrig.res - primate rate4site raw output (this script z score normalizes to match .res output)
		rodent_SE_out - rodent shannon entropy output
		rodent_r4sOrig.res - rodent rate4site raw output (this script z score normalizes to match .res output)
		apn_SE_out - APN shannon entropy output
		apn_r4s.res - APN rate4site z score normalized output
		spike_r4s.res - viral spike rate4site z score normalized output
		spike_binding - viral spike shannon entropy results of binding sites
		spike_no_bin - viral spike shannon entropy results of non binding sites
		all_species.csv - MEME and SLAC datamonkey output for all species
		all_species_minus.csv - MEME and SLAC datamonkey output for all species minus known SARS hosts
		other.csv - MEME and SLAC datamonkey output for the "other mammals" species group
		other_minus.csv - MEME and SLAC datamonkey output for the "other mammals" species group minus known SARS hosts 
		primate.csv - MEME and SLAC datamonkey output for primates
		apn.csv - MEME and SLAC datamonkey output for the APN analysis
		rodent.csv - MEME and SLAC datamonkey output for rodents
		spike.csv - MEME and SLAC datamonkey output for viral spike
		bat.csv - MEME and SLAC datamonkey output for bats
        
    Sequences
        
        apn.fasta - filtered apn sequences for primate apn analysis
        other mammals minus known hosts.fasta - filtered "other mammal sequences" without known sars hosts
        other mammals.fasta - filtered "other mammal sequences"
        primate.fasta - filtered primate sequences
        rodent.fasta - filtered rodent sequences
        spike.fasta - filtered viral spike sequences
        total species minus known hosts.fasta - filtered all species sequences without known sars hosts
        total species.fasta - filtered all species sequences

    Alignments
    
        apn.mafft - mafft aligned primate apn sequences
        other.mafft - mafft aligned "other mammal" sequences
        other_minus.mafft - mafft aligned "other mammal" sequences without known sars hosts
        primate.mafft - mafft aligned primate sequences withough known sars hosts
        rodent.mafft - mafft aligned rodent sequences
        spike.mafft - mafft aligned viral spike sequences
        total.mafft - mafft aligned all species sequences
        total_minus.mafft - codon aligned all species sequences without known hosts
        codon_aligned_apn - codon aligned primate apn sequences
        codon_aligned_other - codon aligned "other mammal" sequences
        codon_aligned_other_minus - codon aligned "other mammal" sequences without known sars hosts
        codon_aligned_primate - codon aligned primate sequences withough known sars hosts
        codon_aligned_rodent - codon aligned rodent sequences
        codon_aligned_spike - codon aligned viral spike sequences
        codon_aligned_total - codon aligned all species sequences
        codon_aligned_total_minus - codon aligned all species sequences without known hosts
        
########################
#     INSTRUCTIONS     #
########################

See master_script.sh for instructions on how to run the code in this repository. The instructions are located in a box at the top of the script.

These instructions are also copied here below:

###############################################################################
#                                                                             #
# download the repository to your local file system                           #
# this master script is provided as a bash script for fully automated running #
# the only variable you need is the path to your directory                    #
#                                                                             #
# enter the repository and type the command in the line below                 #
# run the command: bash master_script.sh path_to_directory                    #
# to run this master script as-is you only need python version 3.8.5          #                                                                               
#                                                                             #
# this will provide all of the output in the main manuscript                  #
# this master script calls the relevant scripts and data in this repository   #
# running this script as-is will give you the exact results in the manuscript #
# if you want to re-build the entire dataset from scratch,                    #
# the commands and programs used are included as comments with instructions   #
#                                                                             #
###############################################################################
