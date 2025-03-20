function [I,Q,fX] = getCC(Y,K,idisp,noise,varargin)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% GETCC calculates channel capacity for an input cell matrix, Y.
%
%
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 



% Data checking
X=cell(size(Y));
for i=1:length(Y)
    A=Y{i};
    if size(A,1)==1 % Unidimensional case: remove NaN individuals
        A=A(:,~isnan(A)); 
    else  % Multidimensional case: remove individuals with any NaNs
        A=A(:,sum(isnan(A))==0); %
    end
    
    % Add random noise to correct for noncontinuous distributions
    A=A+noise*(rand(size(A))-0.5);
    X{i}=A;   
end

% Define dimension and volume
d=size(X{1},1);
V=pi^(d/2)/gamma(d/2+1);



% Loops on loops.
Nq=length(X);
fX=cell(Nq,Nq);
hRS=nan(1,Nq);        
for s1=1:Nq
    for s2=1:Nq            
        Nsamp=size(X{s2},2);
        if s1==s2 %same well            
            [~,DistK]=annquery(X{s2},X{s1},K+1);
            Dk=DistK(K+1,:)';
            fX{s1,s2}=(K)./(Nsamp*V*Dk.^d);                 

            % Conditional Entropies
%             hRS(s1)=-sum(log2(max(fX{s1,s2},eps)))/length(Dk);
            hRS(s1)=-sum(log2(fX{s1,s2}(fX{s1,s2}>eps)))/nnz((fX{s1,s2}>eps));

        else   %different well            
            if isequal(X{s2},X{s1})
                [~,DistK]=annquery(X{s2},X{s1},K+1);
                Dk=DistK(K+1,:)';                
                fX{s1,s2}=(K)./(Nsamp*V*Dk.^d); 
            else              
                [~,DistK]=annquery(X{s2},X{s1},K);
                Dk=DistK(K,:)';
                fX{s1,s2}=(K)./(Nsamp*V*Dk.^d); 
            end              
        end
    end
end
       
if idisp==0
    options = optimoptions('fmincon','Algorithm','active-set','Display','off');
elseif idisp==1
    options = optimoptions('fmincon','Algorithm','active-set','Display','iter');
end

infoF = @(q) -getInfo(q,fX,hRS);
Aeq=ones(1,Nq);
beq=1;
LB=zeros(Nq,1);
UB=ones(Nq,1);  

if numel(varargin)>0
    Q=varargin{1};
    I=-infoF(Q);  
else    
    Qstart=ones(Nq,1);
    Qstart=Qstart/sum(Qstart);
    qfit=fmincon(infoF,Qstart,[],[],Aeq,beq,LB,UB,[],options);
    Q=qfit;
    I=-infoF(qfit);    
end



  


