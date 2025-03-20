function [I,Q,fX] = getCC (Y,K,idisp,varargin)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% getCC calculates channel capacity for an input cell matrix, Y.
%
%
% Y         Cell array (size = 1 x n) of sets of individual responses to n different inputs.
%           Each cell is of size [r,c] -> r dimensions x c individuals (r is same across
%           all cells of Y) 
% K         Use Kth nearest neighbor (passed to aproxmiate KNN algorithm)
% idisp     (true/false) show verbose output
% varargin  (optional) entropies provided
%
%
% fitI      Extrapolated mutual information for entire set
% info1     Mutual information for full set
% I         Mutual information for all subsets
% Q         "Optimal" weights (of all possible inputs) that maximizes mutual information
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Create new array, X, which holds nan-filtered data 
Nsig=length(Y);
X=cell(Nsig,1);
for i=1:Nsig
    A=Y{i};
    if size(A,1)==1 % Unidimensional case: remove NaN individuals
        A=A(:,~isnan(A));
    else % Multidimensional case: remove individuals with any NaNs
        A=A(:,sum(isnan(A))==0); 
    end   
    A=A+(1e-6)*rand(size(A));        
    X{i}=A;   
end


% Define number of dimensions (e.g. timepoints) and volume of unit sphere of dimension=d
d=size(X{1},1); %dimension
V=pi^(d/2)/gamma(d/2+1);

% (STEP 1) Use K-nearest neighbors to calculate: (d-dimensional) probability densities
%    - probability densities (n x n cell matrix [n= # of inputs])
%    - conditional entropies
%
% [n1xn2] conditional probability matrix calculates the probability of getting a response 
% (a query, pulled from response set X{n1}), given another set of reference responses (from X{n2})



Nq=length(X); % # of inputs
fX=cell(Nq,Nq); % Conditional probabilities
hRS=nan(1,Nq); % Conditional entropies       

for s1=1:Nq % Outer loop: "query points" for KNN algorithm; rows of hRS
    for s2=1:Nq  % Inner loop: reference points; columns of hRS          
        Nsamp=size(X{s2},2); 
        % Diagonal of fX (same-well comparisons) -> Probability densities + conditional entropies
        if s1==s2           
            [~,DistK]=annquery(X{s2},X{s1},K+1); % Calculate K-nearest distances for all points
            Dk=DistK(K+1,:)'; % Get distance to (k+1)th nearest neighbor for all samples
            fX{s1,s2}=(K)./(Nsamp*V*Dk.^d); % Formula (2.10) for probability distribution             
            
            if nargin>3 % (optional probabilities are provided)
                if numel(varargin{1})>0  % (ensure entropies are actually a vector)
                    Ftrue=varargin{1};
                    hRS(s1)=-sum(log2(Ftrue{s1,s2}(Ftrue{s1,s2}>eps)))/nnz((Ftrue{s1,s2}>eps));
                else 
                    hRS(s1)=-sum(log2(fX{s1,s2}(fX{s1,s2}>eps)))/nnz((fX{s1,s2}>eps));
                end
            else % Calculate conditional entropies using using eqn 2.12
                hRS(s1)=-sum(log2(fX{s1,s2}(fX{s1,s2}>eps)))/nnz((fX{s1,s2}>eps));
            end
        
        % Off-diagonal (different-well comparision) -> Probability densities only.
        else    
            if isequal(X{s2},X{s1})
                [~,DistK]=annquery(X{s2},X{s1},K+1);
                Dk=DistK(K+1,:)';                
                fX{s1,s2}=(K)./(Nsamp*V*Dk.^d); 
            else              
                [~,DistK]=annquery(X{s2},X{s1},K);
                Dk=DistK(K,:)';
                fX{s1,s2}=(K)./((Nsamp+1)*V*Dk.^d); 
            end              
        end
    end
end
% If probabilities are provided, assign into fX
if nargin>3 
    if numel(varargin{1})>0    
        fX=Ftrue;    
    end
end


% (STEP 2) Optimize input probabilities for maximum channel capacity
% Channel capacities are calculated as the difference between the sum of conditional entropies H(R|S),
% and the total entropy of the response, H(R).
%
%     - H(R|S) = sum(hRS);
%     - H(R) = doubly-weighted sum (see getInfo.m) 


% Set options (toggle verbose/non-verbose display)
if idisp==0
    options = optimoptions('fmincon','Algorithm','active-set','Display','off');
elseif idisp==1
    options = optimoptions('fmincon','Algorithm','active-set','Display','iter');
end

% Define information calculation function (getInfo.m), initialize bounds for optimization.
infoF = @(q) -getInfo(q,fX,hRS);
Aeq=ones(1,Nq);
beq=1;
LB=zeros(Nq,1);
UB=ones(Nq,1);  

if nargin>4 % q's are provided; skip optimization and output.
    if numel(varargin{2})>0
        Qstart=varargin{2};
        I=-infoF(Qstart);  
    end
else % Optimize q for information content, given probability distributions 
    Qstart=ones(Nq,1);
    Qstart=Qstart/sum(Qstart);    
    qfit=fmincon(infoF,Qstart,[],[],Aeq,beq,LB,UB,[],options);
    Q=qfit;
    I=-infoF(qfit);     
end



  


