% osc_feature

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% sin func

% 
% Fs = 12;            % Sampling frequency
% T = 1/Fs;             % Sampling period
% L = 97;             % Length of signal
% t = (0:L-1)*T;        % Time vector % t = 0:1/Fs:8-1/Fs;
% 
% leg_ornot = 0;
% eps = 0.5;
% 
% T = [0.8,1,1.5,2];% 0.5;
% clear x
% for i_x=1:4
% 
%     x(i_x +4,:) = 2+ sin(2*pi/T(i_x)*t);% .*(8-t)./8;
%     % x(i_x +4,:) = x(i_x +4,end:-1:1);
%     x(i_x,:) = x(i_x+4,:)+eps*randn(size(t));
% end
% 
% %
% % x(9,:) = sin(2*pi/T*t)+eps*randn(size(t));
% % x(10,:) = sin(2*pi/T*t);
% time_series = x;
% FramesPerHour= 12;
% 
% sig_stats =get_sig_stats_new(time_series, 'FS', FramesPerHour);
% 
% 
% figure(1)
% for i_x = 1:4
% 
%     subplot(2,2,i_x)
% 
%     plot(t, x(i_x,:),'k','linewidth',1.5);hold on
% 
%     plot(t, x(4+i_x,:),'r','linewidth',1.5);hold on
%     if leg_ornot
%     legend(strcat('op=',num2str(sig_stats.oscpower(i_x),'%.4f'),...
%         ';p2p=',num2str(sig_stats.peak2rms(i_x),'%.4f')),...
%         strcat('op=',num2str(sig_stats.oscpower(4+i_x),'%.4f'),...
%         ';p2p=',num2str(sig_stats.peak2rms(4+i_x),'%.4f')))
%     end
% 
%     title(strcat('T=',num2str(T(i_x)),'hours'))
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% representative cells

load('codon_osc_cells.mat')
time_series = data_tosave.exp{1};
FramesPerHour= 12;

sig_stats =get_sig_stats_v202203(time_series, 'FS', FramesPerHour);


figure(2)
for i_x = 1:12

    subplot(3,4,i_x)

    plot((1:size(time_series,2))/12, time_series(i_x,:),'k','linewidth',1.5);hold on

    plot((1:size(time_series,2))/12, time_series(12+i_x,:),'r','linewidth',1.5);hold on

    legend(strcat('new-op=',num2str(sig_stats.oscpower(i_x),'%.6f'),...
        '; old-op=',num2str(sig_stats.oscpower_old(i_x),'%.4f')...
        ),...
        strcat('new-op=',num2str(sig_stats.oscpower(12+i_x),'%.6f'),...
        '; old-op=',num2str(sig_stats.oscpower_old(12+i_x),'%.4f')...
        ))

    % title(strcat('T=',num2str(T(i_x)),'hours'))
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% fft cal for representative cell


load('codon_osc_cells.mat')
time_series = data_tosave.exp{1};
FramesPerHour= 12;

Fs = 12;            % Sampling frequency
T = 1/Fs;             % Sampling period
L = 97;             % Length of signal
t = (0:L-1)*T;        % Time vector

for i_sig =1:12

    S = time_series(i_sig+12,:);


    X = time_series(i_sig,:);%S + randn(size(t));

    % figure(1)
    % plot(t,X); hold on
    % title('Signal Corrupted with Zero-Mean Random Noise')
    % xlabel('t (hrs)')
    % ylabel('X(t)')

    Y = fft(X);

    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    

    f = Fs*(0:(L/2))/L;

    oscpower = sum(P1(f>=0 & f<=1).^2);
    % figure(2)
    % title('Single-Sided Amplitude Spectrum of X(t)')
    % xlabel('f (Hz)')
    % ylabel('|P1(f)|')

    Y = fft(S);
    P2s = abs(Y/L);
    P1s = P2s(1:L/2+1);
    P1s(2:end-1) = 2*P1s(2:end-1);
oscpowers = sum(P1s(f>=0 & f<=1).^2);

    
    figure(3)
    subplot(3,4,i_sig)
    plot(f,P1,'k','LineWidth',1.5); hold on
    plot(f,P1s,'r','LineWidth',1.5); hold on
    plot(sig_stats.fq,sig_stats.psd(i_sig,:)/10,'k:','LineWidth',2.5); hold on
    plot(sig_stats.fq,sig_stats.psd(i_sig+12,:)/10,'r:','LineWidth',2.5); hold on
    % title('Single-Sided Amplitude Spectrum')
    
%         legend(strcat('osc = ',num2str(oscpower,'%.4f')...
%         ),...
%         strcat('osc = ',num2str(oscpowers,'%.4f')...
%         ))
    xlabel('f (hour^{-1})')
    ylabel('|P1(f)|')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% fft cal for sin

% 
% for i_sig =1:4
% 
%     S = time_series(i_sig+4,:);
% 
%     X = time_series(i_sig,:);%S + randn(size(t));
% 
% 
%     Y = fft(X);
% 
%     P2 = abs(Y/L);
%     P1 = P2(1:L/2+1);
%     P1(2:end-1) = 2*P1(2:end-1);
% 
%     f = Fs*(0:(L/2))/L;
% 
%     % figure(2)
%     % title('Single-Sided Amplitude Spectrum of X(t)')
%     % xlabel('f (Hz)')
%     % ylabel('|P1(f)|')
% 
%     Y = fft(S);
%     P2s = abs(Y/L);
%     P1s = P2s(1:L/2+1);
%     P1s(2:end-1) = 2*P1s(2:end-1);
% 
%  
%     figure(3)
%     subplot(2,2,i_sig)
%     plot(f,P1,'k','LineWidth',1.5); hold on
%     plot(f,P1s,'r','LineWidth',1.5); hold on
%     % title('Single-Sided Amplitude Spectrum')
%     xlabel('f (hour^{-1})')
%     ylabel('|P1(f)|')
%     title(strcat('T=',num2str(T(i_sig)),'hours'))
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fay's osc feature

% addpath('../NFkB_codon/')
% 
% 
% load('codon_osc_cells.mat')
% time_series = data_tosave.exp{1};
% FramesPerHour= 12;
% 
% Fs = 12;            % Sampling frequency
% T = 1/Fs;             % Sampling period
% L = 97;             % Length of signal
% t = (0:L-1)*T;        % Time vector
% 
% peak_range_min =50;
% peak_range_max =140;
% 
% noise_time_diff=35;% if multiple peaks are found within 35 min, most likely due to noise
% noise_peak_diff=0.0088; %An %0.15
% 
% for i_sig = 5%  1:12
%     
%     
%     sig = time_series(i_sig,:);
%     smooth_method = 'moving';
%     smooth_window = 0.03;
%     
%     %first smoothing pass
%     select_cell_smooth = smooth(sig);
%     %another smoothing pass
%     select_cell_smooth_moving_col = smooth(select_cell_smooth, smooth_window, smooth_method);
%     %convert from column to row
%     select_cell_smooth_moving = select_cell_smooth_moving_col.';
%     
%     X  = select_cell_smooth_moving ;
%     
%     
%     sig = time_series(i_sig+12,:);
%     
%     %first smoothing pass
%     select_cell_smooth = smooth(sig);
%     %another smoothing pass
%     select_cell_smooth_moving_col = smooth(select_cell_smooth, smooth_window, smooth_method);
%     %convert from column to row
%     select_cell_smooth_moving = select_cell_smooth_moving_col.';
%     
%     S= select_cell_smooth_moving ;
%     
%     osc_info_S = get_osc_content_trough_peak( S,peak_range_min, peak_range_max,noise_time_diff,noise_peak_diff);
%     osc_info_X = get_osc_content_trough_peak( X,peak_range_min, peak_range_max,noise_time_diff,noise_peak_diff);
%     
%     
%     osc_info_S.sum_first_half_pk_trough_ratio
%     osc_info_X.sum_first_half_pk_trough_ratio
%     
%     osc_info_S.sum_second_half_pk_trough_ratio
%     osc_info_X.sum_second_half_pk_trough_ratio
%     
%     
%     subplot(3,4,i_sig)
%     
%     %     plot((1:size(time_series,2))/12, time_series(i_sig,:),'k','linewidth',1.5);hold on
%     %
%     %     plot((1:size(time_series,2))/12, time_series(12+i_sig,:),'r','linewidth',1.5);hold on
%     
%     plot((1:size(time_series,2))/12, X,'k','linewidth',1.5);hold on
%     
%     plot((1:size(time_series,2))/12, S,'r','linewidth',1.5);hold on
%     
%     legend(strcat('earlyp2t=',num2str(osc_info_X.sum_first_half_pk_trough_ratio),...
%         '; latep2t=',num2str(osc_info_X.sum_second_half_pk_trough_ratio)...
%         ),...
%         strcat('earlyp2t=',num2str(osc_info_S.sum_first_half_pk_trough_ratio),...
%         '; latep2t=',num2str(osc_info_S.sum_second_half_pk_trough_ratio)...
%         ))
% end