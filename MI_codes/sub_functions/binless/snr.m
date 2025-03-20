function [s2n,signal,Noise] = snr(Z)


Nsig = length(Z);
meanZ=zeros(size(Z{1},1),Nsig);
noise=zeros(Nsig,1);
for i=1:Nsig    
    meanZ(:,i)=mean(Z{i},2);
    [~,DistZ]=annquery(meanZ(:,i),Z{i},1);
    noise(i)=sqrt(sum(DistZ.^2));
%     noise(i)=mean(DistZ);    
end

mmeanZ=mean(meanZ,2);
[~,DistmeanZ]=annquery(mmeanZ,meanZ,1);
signal=sqrt(sum(DistmeanZ.^2));
% signal=mean(DistmeanZ);

Noise=mean(noise);
s2n=signal/Noise;
    
    