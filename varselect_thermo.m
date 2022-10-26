
function [binaries, energies, state_out] = varselect_thermo(inputs, data, phis, a, b, atau, btau, way3, target, gibdraws, thermdraws, state_in)

% uses a thermodynamic 'spin flip' Monte Carlo sampler for variable
% selection in a large input dataset

% 'm' is the length of the dataset and 'n' is the number of inputs
[~, n] = size(inputs);

% one basis function for every main effect and two way interaction in this
% version. more terms possibly added later . . . 
vecsize = n + n*(n-1)/2;
if way3
    vecsize = vecsize + n*(n-1)*(n-2)/6;
end

binaries = zeros(thermdraws, vecsize);
energies = zeros(thermdraws, 1);

% get the initial distribution of spins
if nargin == 13
    spinvec = state_in;
else
    spinvec = round(rand(1, vecsize));
end

% the full-sized interaction matrix
mtx = zeros(vecsize, n);
mtx(1:n,1:n) = eye(n);
ind = n+1;
for i=1:n-1
    for j=i+1:n
        mtx(ind, i) = 1;
        mtx(ind, j) = 1;
        ind = ind + 1;
    end
end
if way3
    ind = n + n*(n-1)/2 + 1;
    for i=1:n-2
        for j=i+1:n-1
            for k=j+1:n
                mtx(ind, i) = 1;
                mtx(ind, j) = 1;
                mtx(ind, k) = 1;
                ind = ind + 1;
            end
        end
    end
end

datsize = length(data);
X = zeros(datsize, vecsize);

for i=1:datsize
    
    phind = ceil(inputs(i,:)*499);
    phind = phind + (phind == 0);
    
    for j=1:vecsize
        
        phi = 1;
        for k=1:n
            
            num = mtx(j,k);
            
            if num
                x = 499*inputs(i,k) - phind(k) + 1;
                phi = phi*(phis{num}.zero(phind(k)) + phis{num}.one(phind(k))*x + phis{num}.two(phind(k))*x^2 + phis{num}.three(phind(k))*x^3);
            end
        end
        
        X(i,j) = phi;
        
    end
end

% calculate the initial energy & add it to the distribution along with its
% configuration
binaries(1,:) = spinvec;
damtx = zeros(sum(spinvec),n);
Xin = zeros(datsize, sum(spinvec));
ind = 1;
for j=1:vecsize
    if (spinvec(j))
        damtx(ind,:) = mtx(j,:);
        Xin(:,ind) = X(:,j);
        ind = ind + 1;
    end
end
[~, ~, ~, ~, ~, energies(1)] = gibbs(inputs, data, phis, Xin, damtx, a, b, atau, btau, gibdraws);

temp = abs(energies(1)/2);
accept = 1;

for i=2:thermdraws
    
    flipper = randi(vecsize);
    if spinvec(flipper)
        spinvec(flipper) = 0;
    else
        spinvec(flipper) = 1;
    end
    
    damtx = zeros(sum(spinvec),n);
    Xin = zeros(datsize, sum(spinvec));
    ind = 1;
    for j=1:vecsize
        if (spinvec(j))
            damtx(ind,:) = mtx(j,:);
            Xin(:,ind) = X(:,j);
            ind = ind + 1;
        end
    end
    
    [~, ~, ~, ~, ~, energie] = gibbs(inputs, data, phis, Xin, damtx, a, b, atau, btau, gibdraws);
    
    % Metropolis
    if energie < energies(i-1)
        binaries(i,:) = spinvec;
        energies(i) = energie;
        accept = accept + 1;
    elseif rand < exp((energies(i-1) - energie)/temp)
        binaries(i,:) = spinvec;
        energies(i) = energie;
        accept = accept + 1;
    else
        binaries(i,:) = binaries(i-1,:);
        energies(i) = energies(i-1);
        spinvec = binaries(i,:);
    end
    
    if mod(i, 100) == 0

        if accept/100 < target
            temp = 2*temp;
        elseif accept/100 > target
            temp = 0.5*temp;
        end
        disp(temp)
        accept = 0;

%         if accept/100 < 0.8*target
%             temp = 1.5*temp;
%         elseif accept/100 > 1.2*target
%             temp = 0.5*temp;
%         end
%         disp(temp)
%         accept = 0;
    end
            
        
    disp([i energies(i)])
end

state_out = binaries(thermdraws,:);

end