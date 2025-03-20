clear nfkb collect_feature_vects
load('./example_data_format/mutual_info_cal_data_example.mat')
nfkb = nfkb(1:3);

%nfkb_eg = nfkb;
data_save_file_path_1 = '../raw_data2023/';%_fay_parameter/';

load(strcat(data_save_file_path_1,'Ade_all_stim_unstim_codon.mat'));%All_ligand_codon_2023.mat'))% data,collect_feature_vects,metrics));

color_vec = {[1,0,0],[0,1,0],[0,0,1],[0.7,0.7,0.7]};


if 1
    ligand_all = {'TNF','Pam3CSK','CpG', 'LPS','PolyIC'};
    codon_list = {'TotalActivity','Duration','EarlyVsLate','Speed','PeakAmplitude','OscVsNonOsc'};
    
    
    index_vec = {[1,2,3,16],...%TNF
        [10:12,16],... %Pam3CSK
        [4:6,16],... % CpG
        [7:9,16],... % LPS
        [13:15,16]};%PolyIC
    dose_all = {'Low','Medium','High','Unstim'};
    for i_ligand = 1:length(ligand_all)
        index_data = index_vec{i_ligand};
        
        for i_index_data =  1:length(index_data)%[2,3]%
            num_cells = length(collect_feature_vects.Duration{index_data(i_index_data)});
            x_val = -ones(36*num_cells,1);
            y_val = x_val;
            for i_codon =1:length(codon_list)
                for j_codon = 1:length(codon_list)
                    x_start = (i_codon-1)*6*num_cells + (j_codon-1) * num_cells + 1;
                    y_start = (j_codon-1)*6*num_cells + (i_codon-1) * num_cells + 1;
                    
                    if i_codon == j_codon
                        x_val(x_start:x_start+num_cells - 1 ) = 0;
                        y_val(y_start:y_start+num_cells - 1)= 0;
                    else
                        x_val(x_start:x_start+num_cells - 1 ) = collect_feature_vects.(codon_list{i_codon}){index_data(i_index_data)}+(i_codon-1)*1.2;
                        y_val(y_start:y_start+num_cells - 1)= collect_feature_vects.(codon_list{i_codon}){index_data(i_index_data)}+(i_codon-1)*1.2;
                    end
                end
                
                codon_data = collect_feature_vects.(codon_list{i_codon}){index_data(i_index_data)}+(i_codon-1)*1.2;
                pts = linspace(min(min(codon_data),(i_codon-1)*1.2),max(max(codon_data),i_codon*1.2),100);
                
                
                [f_codon(i_index_data,:,i_codon),xi_codon(i_index_data,:,i_codon)] = ksdensity(codon_data,pts,...
                    'Function','pdf','Bandwidth',0.05);%,'Bandwidth',bw
                
            end
            h(i_index_data) = scatter(x_val,y_val,1,color_vec{i_index_data},'filled');hold on
        end
        
        scale_factor = max(f_codon,[],[1 2]);
        
        for i_index_data = 1:length(index_data)
            
            for i_codon = 1:length(codon_list)
                plot(xi_codon(i_index_data,:,i_codon),f_codon(i_index_data,:,i_codon)/scale_factor(i_codon)+(i_codon-1)*1.2,'Color',color_vec{i_index_data},'LineWidth',1); hold on
            end
            
        end
        
        
        for i_x = 0:6
            plot([i_x*1.2-0.1,i_x*1.2-0.1],[-0.2,6*1.2],'Color',[0.7,0.7,0.7]);hold on
            plot([-0.2,6*1.2],[i_x*1.2-0.1,i_x*1.2-0.1],'Color',[0.7,0.7,0.7]);hold on
            
        end
        
        for i_codon = 1:length(codon_list)
            
            text(i_codon*1.2-0.7,-0.5,codon_list{i_codon},'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize',7)
            
            %             if mod(i_codon,2)
            %            text(i_codon*1.2-0.7,-0.5,codon_list{i_codon},'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize',7)
            %             else
            %               text(i_codon*1.2-0.7,-0.9,codon_list{i_codon},'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize',7)
            %
            %             end
            text(-0.8,i_codon*1.2-0.6,codon_list{i_codon},'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize',7)
            
        end
        text(3,7.5,ligand_all{i_ligand},'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize',9)
        xlim([-2,8])
        ylim([-2,8])
        legend([h(1),h(2),h(3),h(4)],dose_all,'Location','northeast','box','off');
        
        
        
        % medium size figure
        paperpos = [0,0,100,100]*5;
        papersize = [100,100]*5;
        draw_pos = [20,20,60,60]*5;
        set(gcf, 'PaperUnits','points')
        set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)
        %to delete
        fig_save_path = '../SubFigures2023/';
        
        set(gca, 'fontsize',7,'XTick',{},'YTick',{})
        saveas(gcf,strcat(fig_save_path,'codon_distrib_ade_exp_',ligand_all{i_ligand}),'epsc')
        close()
    end
end