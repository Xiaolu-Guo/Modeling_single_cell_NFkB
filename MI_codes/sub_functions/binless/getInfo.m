function [f,Hr] = getInfo(q,Fcond,hRS)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% GETINFO calculates mutual information for:
% 
% q         a set of input probabilities (sum(q)=1)
% Fcond     an [nxn] conditional probability matrix
% hRS       a vector of conditional entropies
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 
% Get rows of conditional probability matrix - one "query", all "references"
% Sum across this row, individually weighting each point by appropriate
% weight.
Nq=length(Fcond);
F=cell(Nq,1); 
Wq=cell(Nq,1);
for k=1:Nq
    F{k}=cat(2,Fcond{k,:})*q; 
    Wq{k}=ones(size(F{k}))*q(k)/sum(F{k}>eps);              
end

% Combine all probabilities/weights
F_all=cat(1,F{:});
Wq=cat(1,Wq{:});

%F_all(F_all<=eps)=nan; 
% Unweight any zero-values probability estimates
Wq(F_all<=eps) = 0;

% Calculate non-conditional entropy (via eqn 2.13)
Hr = -sum(log2(F_all).*Wq); 

% Entropy = difference of H(R) and H(R|S), summed over the probabilities of each response. 
f=Hr-hRS*q;    
    