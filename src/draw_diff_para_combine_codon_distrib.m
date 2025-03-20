function [] = draw_diff_para_combine_codon_distrib(data_save_file_path,fig_save_path)
%This is for visualizing Monolix estimation
% all; -66,-99,-100,-101,-99-100-101,-receptor1,-rpt2,-rcpt3
data_proj_nums ={[18,22,21,20,19,16,24,23];%TNF
    [25,30,34,35,36,29,31,33,32];%LPS
    [26,41,42,43,44,38,40,39];% CpG
    [28,46,47,48,49,45,51,50];
    [27,55,56,57,58,59,52,53,54]};% P3CSK

ligand_all = {'TNF';'LPS';'CpG';'Pam3CSK';'polyIC'};%polyIC

fig_opt.paper_opt.paperpos=[0,0,150,180]*3;
fig_opt.paper_opt.papersize=[150 180]*3;

for i_ligand = 1:length(data_proj_nums)
    data_proj_num_vec = data_proj_nums{i_ligand};    

    load(strcat(data_save_file_path,ligand_all{i_ligand},'_diff_para_combine_codon.mat'))% data,collect_feature_vects,metrics));
       
        %% codon distribution
    % violin_plot_codon_seperate(collect_feature_vects,fig_save_path) % draw each
    % codon and sti seperately, and save
    fig_opt.save_file = strcat(fig_save_path,ligand_all{i_ligand},'_diff_para_combine_codon_distrib');
    violin_plot_codon(collect_feature_vects,fig_opt) % draw all codon and sti in one

 
end



