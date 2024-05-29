
function delta = bss_derivatives(x, betas, phis, mtx, derv, range)

% x are inputs
% betas are coefficients

% phis are the spline coefficients for the basis functions (cell array)
 
% mtx is the 'interaction matrix' -- a matrix each row of which corresponds
% to a term in the expansion, each column corresponds to an input. if the
% column is zero there's no corresponding basis function in the term; if
% it's greater than zero it corresponds to the order of the basis function

% derv is the order of differentiation -- 1 for a first derivative, 2 for a
% second. it's a vector the same size as the column dimensions of 'mtx' and
% the numbers go in the place of the input that should be differentiated
% with respect to

% range is the range of the input data used in the normalization, such that
% derivatives can be converted from the normalized scale to the actual
% scale. it's a vector the same size as the column dimensions of 'mtx'
    
[m,n] = size(mtx);

if length(derv) ~= n
    ME = MException('bss_derivatives:wrongsize', 'derv must be equal in size to the number of inputs');
    throw(ME)
end

for i=1:n
    if derv(i) ~= 0 && derv(i) ~= 1 && derv(i) ~= 2
        ME2 = MException('bss_derivatives:wrongorder', 'only valid differentiation orders are 1 and 2');
        throw(ME2)
    end
end

if sum(derv) == 0
    error('bss_derivatives: must have at least one derivative. for function eval without derivatives use ''bss_eval''');
end

[mx, ~] = size(x);
[mbet, ~] = size(betas);

delta = zeros(mx,mbet);

if isstruct(phis{1})
    bsplines = true;
else
    bsplines = false;
    phind = ceil(x*499);

    set = (phind == 0);
    phind = phind + set;

    r = 1/499;
    xmin = (phind-1)*r;
    X = (x-xmin)/r;
end

for i=1:m

    phi = ones(mx,1);
    for j=1:n

        num = mtx(i,j);

        if num
            derp = derv(j);
            if derp == 0
                if bsplines
                    phi = phi.*fnval(phis{num},x(:,j));
                else
                    phi = phi.*(phis{num}.zero(phind(:,j)) + phis{num}.one(phind(:,j)).*X(:,j) + ...
                        phis{num}.two(phind(:,j))*X(ii,j).^2 + phis{num}.three(phind(:,j)).*X(:,j).^3);
                end
            elseif derp == 1
                if bsplines
                    phi = phi.*fnval(phis{num+500},x(:,j))/range(j);
                else
                    phi = phi.*(phis{num}.one(phind(:,j)) + 2*phis{num}.two(phind(:,j)).*X(:,j) + ...
                        3*phis{num}.three(phind(:,j)).*X(:,j).^2)./(range(j)/499);
                end
            elseif derp == 2
                if bsplines
                    phi = phi.*fnval(phis{num+1000},x(:,j))/range(j)^2;
                else
                    phi = phi.*(2*phis{num}.two(phind(:,j)) + 6*phis{num}.three(phind(:,j)).*X(ii,j))/(range(j)/499)^2;
                end
            end
        elseif derv(j)
            phi = zeros(mx,1);
            break;
        end

    end

    delta = delta + repmat(phi,1,mbet).*repmat(betas(:,i+1)',mx,1);

end






end