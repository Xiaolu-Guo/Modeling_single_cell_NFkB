%plot distribution
% 
filepath='./'

estimates=readtable(strcat(filepath,'Model_all_data_core_parameter_results_0423.txt'));

color_num=[1,1,1,1,2,1,2,2,2,2,2,2,2,2];
color_num=[1,1,1,1,1,1,1,1,1,1,1,1,1,1];
color_name={'r','b'};
%range=[0,1];

%x_mean=0;

for ii=1:11
    range=[estimates.min(ii),estimates.max(ii)];

x_mean=estimates.mean(ii);

y_std=estimates.std(ii);

mu=log((x_mean-range(1))./(range(2)-x_mean));

%mu=x_mean;

sigma=y_std;


x_vec=linspace(range(1),range(2),2000);


x=x_vec(2:end-1);

pdf_cal=1/(sigma*sqrt(2*pi))*(range(2)-range(1))./(x-range(1))./(range(2)-x).*exp(-(log((x-range(1))./(range(2)-x))-mu).^2./(2*sigma^2));

y=log((x-range(1))./(range(2)-x));

pd=makedist('Normal','mu',mu,'sigma',sigma);

pdf_logit_normal=pdf(pd,y);
cdf_logit_normal=cdf(pd,y);

figure(1)
subplot(6,2,ii)
plot(x,(range(2)-range(1))./((x-range(1)).*(range(2)-x)).*pdf_logit_normal,'LineWidth',2,'Color',color_name{color_num(ii)});hold on
%plot(x,pdf_cal,'LineWidth',2);hold on
title(estimates.params{ii})
yl=ylim();
plot([x_mean,x_mean],yl,'--','LineWidth',2,'Color',[0.5,0.5,0.5]); hold on

end


ii=12;
 
x_mean=estimates.mean(ii);

y_std=estimates.std(ii);

mu=x_mean;

sigma=y_std;

x_vec=linspace(mu-3*sigma,mu+3*sigma,2000);
x=x_vec(2:end-1);
pd=makedist('Normal','mu',mu,'sigma',sigma);

pdf_normal=pdf(pd,x);
cdf_normal=cdf(pd,x);

figure(1)
subplot(6,2,ii)
plot(x,pdf_normal,'LineWidth',2,'Color',color_name{color_num(ii)});hold on
%plot(x,pdf_cal,'LineWidth',2);hold on
title(estimates.params{ii})

