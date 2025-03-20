
names = {'NFkBn','IkBaNFkBn','IKK'};%
% names = {'TAK1'};%
paper_pos = [0,0,100,80]*1.5;
paper_size = [100,80]*1.5;
font_size = 7;

% names = {'TAK1'};%
% names = {'IKK'};%

sti_vec = {'TNF','LPS','Pam3CSK','CpG','polyIC'};
%% set the ligand, dose val, dose label , and dose scale.
for i_sti = 1:length(sti_vec)
    sti = sti_vec{i_sti};
    switch sti
        case 'TNF' % ligand
            doses = [0.33,3.3,33];% dose val
            dose_str = {'0.33ng/mL','3.3ng/mL','33ng/mL'};% dose label
            dose_scale = 1/5200.0; % dose scale, DO NOT CHANGE for any ligand. This is the scale factor to change the ligand into uM.
        case 'LPS'
            doses = [0.33,3.3,33];
            dose_str = {'0.33ng/mL','3.3ng/mL','33ng/mL'};
            dose_scale = 1/24000.0;
        case 'Pam3CSK'
            doses = [10,33,100];
            dose_str = {'10ng/mL','33ng/mL','100ng/mL'};
            dose_scale = 1/1500.0;
        case 'CpG'
            doses = [33,100,330];
            dose_str = {'33nM','100nM','330nM'};
            dose_scale = 1/1000.0;
        case 'polyIC'
            doses = [3.3e3,33e3,100e3];
            dose_str = {'3.3ug/mL','33ug/mL','100ug/mL'};
            dose_scale = 1/5000000.0;
    end
    
    %% initialization
    options = struct;
    options.DEBUG = 1;
    options.SIM_TIME = 8*60;
    [v0.PARAMS, v0.SPECIES] = nfkbInitialize();%_TNF_TRAF2_paraset2
    options.v.PARAMS = v0.PARAMS;
    options.v.SPECIES = v0.SPECIES;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Simulate all doses (only need to equilibrate on first iteration)
    output = [];
    
    for i = 1:length(doses)
        if isempty(output)
            [t,x,simdata] = nfkbSimulate({sti,doses(i)*dose_scale},names, [], {},options);% _TNF_TRAF2
        else
            options.STEADY_STATE = simdata.STEADY_STATE;
            [t,x] = nfkbSimulate({sti,doses(i)*dose_scale},names, [], {},options); % _TNF_TRAF2
        end
        output = cat(3,output,x);
    end
    
    NFkBn_curves = squeeze(output(:,strcmp(names,'NFkBn'),:));
    IkBaNFkBn_curves = squeeze(output(:,strcmp(names,'IkBaNFkBn'),:));
    NFkB_nuc = NFkBn_curves + IkBaNFkBn_curves;
    % NFkBn_curves = squeeze(output(:,strcmp(names,'TAK1'),:));
    IKK_curves = squeeze(output(:,strcmp(names,'IKK'),:));
    
    
    figure(1)
    set(gcf, 'PaperUnits','points')
    set(gcf, 'PaperPosition', paper_pos,'PaperSize',paper_size )%,'Position',draw_pos [20,20,280,280]
    
    for i_dose = 1:size(NFkB_nuc,2)
        plot([0:480]/60,NFkB_nuc(:,i_dose),'LineWidth',1);hold on
    end
    
    legend(dose_str)
    legend('boxoff')
    set(gca,'FontSize',font_size,'FontName','Arial')
    
    title(sti,'FontWeight','b','FontSize',9,'FontName','Arial')
    ylabel('nuclear NFkB (\mu M)','FontWeight','b','FontSize',9,'FontName','Arial')
    xlabel('Time (hrs)','FontWeight','b','FontSize',9,'FontName','Arial')
    
    xlim([0,8])
    saveas(gcf,strcat(fig_save_path,sti,'_NFkB_representative'),'epsc')
    close
    
    figure(2)
    
    set(gcf, 'PaperUnits','points')
    set(gcf, 'PaperPosition', paper_pos,'PaperSize',paper_size )%,'Position',draw_pos [20,20,280,280]
    
    plot([0:480]/60,IKK_curves,'LineWidth',1)
    
    % legend(dose_str)
    set(gca,'FontSize',font_size,'FontName','Arial')
    
    title(sti,'FontWeight','b','FontSize',9,'FontName','Arial')
    ylabel('pIKK (\mu M)','FontWeight','b','FontSize',9,'FontName','Arial')
    xlabel('Time (hrs)','FontWeight','b','FontSize',9,'FontName','Arial')
    xlim([0,8])
    %draw_opt.save_filename ='./fig';
    saveas(gcf,strcat(fig_save_path,sti,'_IKK_representative'),'epsc')
    close
end