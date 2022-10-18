
function [betas, sigs, taus, betahat, X, ev] = gibbs(inputs, data, phis, Xin, discmtx, a, b, atau, btau, draws)

% this version of the sampler increases efficiency by accepting an set of
% inputs 'Xin' derived from earlier iterations. 

% 'inputs' is the set of normalized inputs -- both parameters and model inputs -- with
% columns corresponding to inputs and rows the different
% experimental designs

% 'data' are the experimental results: column vector, with entries
% corresponding to rows of 'inputs'

% 'phis' are a data structure with the spline coefficients for the basis
% functions, built with 'spline_coefficients.txt' and 'splineconvert'

% 'discmtx' is the interaction matrix for the bss-anova function -- rows
% are terms in the function and columns are inputs (cols should line up
% with cols in 'inputs'

% 'a' and 'b' are the parameters of the ig distribution for the
% observation error variance of the data

% 'atau' and 'btau' are the parameters of the ig distribution for the 'tau
% squared' parameter: the variance of the beta priors

% 'draws' is the total number of draws


%% first build the matrix by calculating the corresponding basis function outputs for each set of inputs
[minp, ninp] = size(inputs);

[mmtx, ~] = size(discmtx);
[~, nxin] = size(Xin);

X = [Xin zeros(minp, mmtx-nxin) ones(minp,1)]; % number of data points by number of terms in the function
%X(:,mmtx+1) = ones(minp,1); % representing the beta-naught term

for i=1:minp
    
    phind = ceil(inputs(i,:)*499);
    phind = phind + (phind == 0);
    
    for j=nxin+1:mmtx
        
        phi = 1;
        for k=1:ninp
            
            num = discmtx(j,k);
            
            if num
                x = 499*inputs(i,k) - phind(k) + 1;
                phi = phi*(phis{num}.zero(phind(k)) + phis{num}.one(phind(k))*x + phis{num}.two(phind(k))*x^2 + phis{num}.three(phind(k))*x^3);
            end
        end
        
        X(i,j) = phi;
        
    end
end

% initialize sigsqd and tausqd at the mode of their priors
sigsqd = b/(1+a);
tausqd = btau/(1+atau);

XtX = X'*X;

Xty = X'*data;

[Q, Lamb] = eig(XtX);

Lamb_inv = diag(1./diag(Lamb));

betahat = Q*Lamb_inv*Q'*Xty;
squerr = norm(data - X*betahat)^2;

astar = a + length(data)/2 + (mmtx + 1)/2;
atau_star = atau + (mmtx+1)/2;

dtd = data'*data;

%% Gibbs iterations

betas = zeros(draws, mmtx+1);
sigs = zeros(draws, 1);
taus = zeros(draws, 1);

lik = zeros(draws, 1);
n = length(data);

for k = 1:draws
    
    Lamb_tausqd = Lamb + (1/tausqd)*eye(mmtx+1);
    Lamb_tausqd_inv = diag(1./diag(Lamb_tausqd));

    mun = Q*Lamb_tausqd_inv*Q'*Xty;
    S = Q*diag(sqrt(diag(Lamb_tausqd_inv)));

    vec = normrnd(0,1,[mmtx+1,1]);
    betas(k,:) = (mun + sqrt(sigsqd)*S*vec)';

    lik(k) = -(n/2)*log(sigsqd) - (squerr + (betahat' - betas(k,:))*XtX*(betahat - betas(k,:)'))/(2*sigsqd);
    %lik(k) = -(n/2)*log(sigsqd) - norm(data - X*betas(k,:)')^2/(2*sigsqd);
    %evs(k) = sqrt(prod(diag(Lamb_tausqd)))/(2*pi*sigsqd)^(mmtx/2 + 0.5);

    vecc = mun - betas(k,:)';    
    bstar = b + 0.5*vecc'*((XtX + (1/tausqd)*eye(mmtx+1))*vecc) + 0.5*dtd - 0.5*mun'*Xty;
    sigsqd = 1/gamrnd(astar, 1/bstar);
    
    sigs(k) = sigsqd;
    btau_star = (1/(2*sigsqd))*(betas(k,:)*betas(k,:)') + btau;
    tausqd = 1/gamrnd(atau_star, 1/btau_star);
    taus(k) = tausqd;
    
end

%% calculate the evidence

ev = (mmtx+1)*log(n) - 2*max(lik);

X = X(:,1:mmtx);

end