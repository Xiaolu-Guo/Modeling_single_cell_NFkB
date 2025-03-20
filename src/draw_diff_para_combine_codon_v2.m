function [] = draw_diff_para_combine_codon_v2(data_save_file_path,fig_save_path)
%This is for visualizing Monolix estimation
% all; -66,-99,-100,-101,-99-100-101,-receptor1,-rpt2,-rcpt3
data_proj_nums ={[18,22,21,20,19,16,24,23];%TNF
    [25,30,34,35,36,29,31,33,32];%LPS
    [26,41,42,43,44,38,40,39];% CpG
    [28,46,47,48,49,45,51,50];
    [27,55,56,57,58,59,52,53,54]};% P3CSK
% data_proj_nums = {[18,16,24];% TNF
%     [25,29,32]};% LPS


ligand_all = {'TNF';'LPS';'CpG';'Pam3CSK';'polyIC'};%polyIC


paperpos=[0,0,35,90]*3;
papersize=[35 90]*3;
draw_pos=[10,10,30,80]*3;

ligand_index = [1;2];
data_proj_index = {[1,6,7];
    [1,6,9]};

for i_ligand = 1:length(ligand_index)
    data_proj_num_vec = data_proj_nums{i_ligand};    

    load(strcat(data_save_file_path,ligand_all{ligand_index(i_ligand)},'_diff_para_combine_codon.mat'))% data,collect_feature_vects,metrics));
       
    %% codon comparison between different para combination
    figure(i_ligand)
    
    set(gcf, 'PaperUnits','points')
    set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)
    
    % compare codons
    codon_name = {'OscVsNonOsc','PeakAmplitude','EarlyVsLate','TotalActivity','Speed','Duration'};
    codon_diff = struct();
    
    for i_codon = 1:6
        codon_diff.(codon_name{i_codon}) = zeros(length(data_proj_num_vec),1);
        codon_diff_v0.(codon_name{i_codon}) = zeros(length(data_proj_num_vec),max(cell2mat(data.info_dose_index)));
        kk=1;
        
        for ii=1:length(data_proj_num_vec)
            codon_diff_val = [];
            for jj=1:max(cell2mat(data.info_dose_index))
                codon_diff_val = [codon_diff_val;abs(collect_feature_vects.(codon_name{i_codon}){2*kk-1}-collect_feature_vects.(codon_name{i_codon}){2*kk})];
                codon_diff_v0.(codon_name{i_codon})(ii,jj) = mean( abs(collect_feature_vects.(codon_name{i_codon}){2*kk-1}-collect_feature_vects.(codon_name{i_codon}){2*kk}));
                kk=kk+1;
            end
            codon_diff.(codon_name{i_codon})(ii) = mean(codon_diff_val);
            %codon_diff.(codon_name{i_codon})(ii) = mean(codon_diff_v0.(codon_name{i_codon})(ii,:));

        end
        % codon_all = sum(codon_diff.(codon_name{i_codon}),2);
        % codon_all <= codon_all(2)
        
        subplot(6,1,i_codon)
        
        bar(codon_diff.(codon_name{i_codon})(data_proj_index{i_ligand}),'stacked')
        % ylabel(codon_name{i_codon})
        set(gca, 'XTick',[], 'XTickLabel',{})
        
        set(gca,'fontsize',8,'fontweight','b')
        set(gca, 'FontName', 'Arial')
    end

    saveas(gcf,strcat(fig_save_path,ligand_all{ligand_index(i_ligand)},'_diff_para_combine_codon_comparison_delay_receptor_deg'),'epsc')
    close
 
end



