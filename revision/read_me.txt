For machine learning applied on all ligand doses:
1. Run transfer_data_ML_alldoses_exp_sim_20241025.m, to save the experimental data and coding data in the machine learning format, under ./data_signaling_codons_Machine_learning_format/. 
2. Run data_signaling_codons_Machine_learning_format/cross_exp_fitting.ipynb on jupyter notebook, with python 3 (vers 3.12.4); this will apply random forest classification, with 5 fold cross validation, with parameter grid search optimized. confus_mat will be saved under ./data_signaling_codons_Machine_learning_format/. 
3. Run ML_15_Stim_exp_fitting.m to plot the matrix, using confusion_mat saved from last step, and save it to ../../SubFigures2023/

For NLME fitting process:
1. Downsampling the dataset: 
1-1) Downsampling Ade's data, for Pam and pIC condition, run monolix
1-2) todo: and check the fitting performance. 
Calculate the signaling codon, using all 19(?) conditions, check W-dist, RMSD ratio, and classification.

2. To do: Rerun all the fitting:
2-1) rescaling the data, using only 3 doses for CpG and LPS; change the cell numbers; generating datafile with only the 15 conditions.
2-2) rewrite the model file, the parameter initial value file. 
2-3) run it. 

For specific response vs confusion response:
1. Run run_me_revision_nc.m, which called the following functions
1-1) draw_single_cell_confusion_revision_20241027.m
This function will generate single-cell signaling codon difference distributions for each signaling codon and for each stimulus pair; 
This function will also compare within single-cell responses, the specific response integral vs confusion response integral, and peak.
scSRS_corr_heterogeneity_202410.m
This function will calculate "prob of perfect SRS vs cell #" curve, under the setting  that only specific cells are more than confusion cells, the group are considered as achieving specificity. 