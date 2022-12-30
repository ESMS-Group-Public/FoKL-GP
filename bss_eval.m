
function delta = bss_eval(x, betas, phis, mtx)

% x are normalized inputs -- rows experiments & columns variables

% betas are coefficients -- columns correspond to terms in the model (rows
% in mtx) and rows to different draws (different models). the last column
% is the beta-naught (constant) term

% phis are the spline coefficients for the basis functions (cell array)
 
% mtx is the 'interaction matrix' -- a matrix each row of which corresponds
% to a term in the expansion, each column corresponds to an input. if the
% column is zero there's no corresponding basis function in the term; if
% it's greater than zero it corresponds to the order of the basis function 

[m,n] = size(mtx);
[mx, ~] = size(x);
[mbet, ~] = size(betas);

delta = zeros(mx,mbet);
phind = ceil(x*499);

set = (phind == 0);
phind = phind + set;

r = 1/499;
xmin = r*(phind-1);
X = (x-xmin)/r;

for ii=1:mx

    for i=1:m

        phi = 1;
        for j=1:n

            num = mtx(i,j);

            if num
                phi = phi*(phis{num}.zero(phind(ii,j)) + phis{num}.one(phind(ii,j))*X(ii,j) + phis{num}.two(phind(ii,j))*X(ii,j)^2 + phis{num}.three(phind(ii,j))*X(ii,j)^3);
            end

        end

        delta(ii,:) = delta(ii,:) + phi*betas(:,i+1)';

    end

    delta(ii,:) = delta(ii,:) + betas(:,1)';

end

end