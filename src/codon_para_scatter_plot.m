% scatter(para_vals(:,2),sig_codons_model(:,1))
figure(3)
x = para_vals(:,5);
y = sig_codons_model(:,2);
% ft = fittype( 'Hill_line_neg(x,n,K,A)' );
% % f = fit( x, y, ft, 'StartPoint', [1, mean(x), 1] );
% f = fit( x, y, ft, 'StartPoint', [1, mean(x), 1] );

ft = fittype( 'Hill_line_pos(x,n,K,A)' );
f = fit( x, y, ft, 'StartPoint', [2, mean(x), 1] );

% f=fit(x,y,'smoothingspline');%  lowess
h = plot(f,x,y);

set(gca,'XScale','log','FontWeight','b','FontSize',14)
% xlabel('TAK1 activiation')%TAK1 activiation 
% xlabel('Time delay')
xlabel('TLR-ligand deg')%TAK1 activiation 
% xlabel('Traffic to Endosome')%TAK1 activiation 
% xlabel('TLR synthesis')%TAK1 activiation 

% ylabel('Oscillation')
ylabel('EarlyVsLate')
% ylabel('Speed')