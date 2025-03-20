function [] = plot_est_distribution_seperate_dose(para_sample, para_name,estimates)
% load distribute_regulator_unimodal_ind_codon.mat

[S,AX,BigAx,H,HAx] = plotmatrix(para_sample(:,:));
ax=cell(size(para_sample,2),size(para_sample,2));

figure(1)
for i=1:size(para_sample,2)
    for j=1:size(para_sample,2)
        
        AX(i,j).XScale='log';
        AX(i,j).YScale='log';
        ax{i,j}.Position = AX(i,j).Position;
        ax{i,j}.Position(3)=ax{i,j}.Position(3)*0.8;
        ax{i,j}.Position(4)=ax{i,j}.Position(4)*0.8;
        
    end
end
hold on

ax_Lim = [0.001, 1];
ax_Tick = 10.^(-10:10);

cell_index = cumsum(estimates.cell_num);
cell_index = [0; cell_index];

figure(2)
for i=1:length(para_name)
    for j=1:length(para_name)
        
        subplot('Position',ax{i,j}.Position)
        
        [~,edges_i] = histcounts(log10(para_sample(:,i)));
        % [~,edges_j] = histcounts(log10(R(:,j)));
        
        if i==j
            
            for i_dose = 1: length(estimates.cell_num)
                histogram(para_sample((cell_index(i_dose)+1):cell_index(i_dose+1),i),10.^edges_i,'Normalization','pdf');hold on
                fitted_para_doses{i_dose} = para_sample((cell_index(i_dose)+1):cell_index(i_dose+1),i);
            end
                    
            H = 0;
            for i_dose = 1:length(fitted_para_doses)
                [h,p] = kstest2(fitted_para_doses{i_dose},fitted_para_doses{mod(i_dose,length(fitted_para_doses))+1})
                H = H +h;
            end
            H
            
            logit_plot(estimates.min(i),estimates.max(i),estimates.mean(i),estimates.std(i));
        else
            
            for i_dose = 1:length(estimates.cell_num)
                scatter(para_sample((cell_index(i_dose)+1):cell_index(i_dose+1),j),para_sample((cell_index(i_dose)+1):cell_index(i_dose+1),i),2,estimates.dose_color{i_dose},'filled');hold on
            end
            
        end
        axc=gca;
        
        XL(2) = prctile(para_sample(:,j),100)*1.3;
        XL(1) = prctile(para_sample(:,j),0)*0.7;
        Xtick = ax_Tick(XL(1)<ax_Tick & ax_Tick<XL(2));
        if length(Xtick)<3
            Xtick = Xtick;
        else
            Xtick =Xtick([1,ceil(end/2),end]);
        end
        
        YL(2) = prctile(para_sample(:,i),100)*1.3;
        YL(1) = prctile(para_sample(:,i),0)*0.7;
        Ytick = ax_Tick(YL(1)<ax_Tick & ax_Tick<YL(2));
        if length(Ytick)<3
            Ytick = Ytick;
        else
            Ytick =Ytick([1,ceil(end/2),end]);
        end
        
        axc.XLim=XL;
        axc.XScale='log';
        
        if i~=j
            
            axc.YLim=YL;
            axc.YScale='log';
            axc.YTick =[];
        end
        
        axc.XTick =[];
        
        if j==1
            
            if i==1
                ylabel('probability','fontweight','b')
            else
                ylabel(para_name{i},'fontweight','b')
                axc.YTick = Ytick ;
            end
        end
        if i==length(para_name)
            xlabel(para_name{j},'fontweight','b')
            axc.XTick = Xtick ;
            
        end
        %
    end
end
hold on

figure(1)
close
end


function [] = logit_plot(estimates_min,estimates_max,estimates_mean,estimates_std)
range=[estimates_min,estimates_max];
x_mean=estimates_mean;
y_std=estimates_std;
mu=log((x_mean-range(1))./(range(2)-x_mean));
%mu=x_mean;

sigma=y_std;
x_vec=linspace(range(1),range(2),2000);
x=x_vec(2:end-1);

try
    pd=makedist('Normal','mu',mu,'sigma',sigma);
    y=log((x-range(1))./(range(2)-x));
    pdf_logit_normal=pdf(pd,y);
    % cdf_logit_normal=cdf(pd,y);
    plot(x,(range(2)-range(1))./((x-range(1)).*(range(2)-x)).*pdf_logit_normal,'LineWidth',2,'Color','k');hold on
catch
    pdf_cal=1/(sigma*sqrt(2*pi))*(range(2)-range(1))./(x-range(1))./(range(2)-x).*exp(-(log((x-range(1))./(range(2)-x))-mu).^2./(2*sigma^2));
    plot(x,pdf_cal,'LineWidth',2,'Color','k');hold on
end

% yl=ylim();
% plot([x_mean,x_mean],yl,'--','LineWidth',2,'Color',[0.5,0.5,0.5]); hold on
end
