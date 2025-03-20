function [responder_ratio_by_peak, responder_ratio] = calculation_responder_5sigma_Supriyadata()
    load('Supriya_data_metrics_2024.mat')

   a = 1; 

       cal_responder_ratio = 1;
    index_single_ligand = {[4];% 1-TNF
        [6,24]; %2-
        [1,10]; %3
        [7,12,15,17,22]; %4
        [2,11,23,27]}; %5
    index_single_all7ligand = {[4];% 1-TNF
        [6,24]; %2-
        [1,10]; %3
        [7,12,15,17,22]; %4
        [2,11,23,27];
        [28]}    
    
    baseline_vals = [];
    for i_index_single_ligand = 1:length(index_single_ligand)
        index_single_vec = index_single_ligand{i_index_single_ligand};
        for i_index_single_vec = 1:length(index_single_vec)
            baseline_vals = [baseline_vals;data.exp{i_index_single_vec}(:,1)];
        end
    end
    
    base_std = nanstd(baseline_vals);
    base_mean = nanmean(baseline_vals);
        
    integral_pos_threshold = (base_mean + base_std*5)*0.5;
    peak_threshold = base_mean + base_std*5;
    
    index_scatter_y = 1;
    
    for i_cond = 1:length(index_single_all7ligand)
        index_exp_ligand = index_single_all7ligand{i_cond};
        responder_ratio(i_cond) = 0;
        responder_ratio_by_peak(i_cond) = 0;
        for i_rpc = 1:length(index_exp_ligand)
            responder_ratio_scatter_y(index_scatter_y) = sum(metrics{index_exp_ligand(i_rpc)}.integrals_pos(:,end) > integral_pos_threshold)/length(metrics{index_exp_ligand(i_rpc)}.integrals_pos(:,end));
            responder_ratio_scatter_x(index_scatter_y) = i_cond;
            responder_ratio(i_cond) =responder_ratio(i_cond)+ responder_ratio_scatter_y(index_scatter_y);
            
            responder_ratio_by_peak_scatter_y(index_scatter_y) = sum(metrics{index_exp_ligand(i_rpc)}.max_amplitude(:,end) > peak_threshold)/length(metrics{index_exp_ligand(i_rpc)}.max_amplitude(:,end));
            responder_ratio_by_peak_scatter_x(index_scatter_y) = i_cond;
            responder_ratio_by_peak(i_cond) = responder_ratio_by_peak(i_cond)+responder_ratio_by_peak_scatter_y(index_scatter_y);
            index_scatter_y = index_scatter_y+1;
        end
        
        if responder_ratio(i_cond)
            responder_ratio(i_cond) = responder_ratio(i_cond)/length(index_exp_ligand);
            responder_ratio_by_peak(i_cond) = responder_ratio_by_peak(i_cond)/length(index_exp_ligand);
        end
    end
    
    responder_ratio_by_peak
    responder_ratio

    
end