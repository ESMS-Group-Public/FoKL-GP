
function [betas, mtx, evs] = emulator(inputs, data, phis, a, b, atau, btau, tolerance, draws, gimmie, way3, aic)

% this version uses the 'Xin' mode of the gibbs sampler

% builds a single-output bss-anova emulator for a stationary dataset in an
% automated fashion

% function inputs:

% 'inputs' is the set of inputs normalized on [0,1]: matrix with
% columns corresponding to inputs and rows the different
% experimental designs

% 'data' are the output dataset used to build the function: column vector, with entries
% corresponding to rows of 'inputs'

% 'phis' are a data structure with the spline coefficients for the basis
% functions, built with 'spline_coefficient.txt' and 'splineconvert' or
% 'spline_coefficient_500.txt' and 'splineconvert500' (the former provides
% 25 basis functions: enough for most things -- while the latter provides
% 500: definitely enough for anything)

% 'a' and 'b' are the shape and scale parameters of the ig distribution for the
% observation error variance of the data. the observation error model is white noise
% choose the mode of the ig distribution to match the noise in the output dataset
% and the mean to broaden it some

% 'atau' and 'btau' are the parameters of the ig distribution for the 'tau
% squared' parameter: the variance of the beta priors is iid normal mean zero with
% variance equal to sigma squared times tau squared. tau squared must be scaled in the
% prior such that the product of tau squared and sigma squared scales with the output
% dataset

% 'tolerance' controls how hard the function builder tries to find a better
% model once adding terms starts to show diminishing returns. a good
% default is 3 -- large datasets could benefit from higher values

% 'draws' is the total number of draws from the posterior for each tested 

% 'gimmie' is a boolean causing the routine to return the most complex
% model tried instead of the model with the optimum bic

% 'way3' is a boolean for turning on or off 3-way interactions

% 'aic' is a boolean specifying the use of the aikaike information
% criterion

% function outputs:

% 'betas' are a draw from the posterior distribution of coefficients: matrix, with
% rows corresponding to draws and columns corresponding to terms in the GP

% 'mtx' is the basis function interaction matrix from the
% best model: matrix, with rows corresponding to terms in the GP (and thus to the 
% columns of 'betas' and columns corresponding to inputs. a given entry in the 
% matrix gives the order of the basis function appearing in a given term in the GP.
% all basis functions indicated on a given row are multiplied together.
% a zero indicates no basis function from a given input is present in a given term

% 'ev' is a vector of BIC values from all of the models
% evaluated

% 

% 'n' is the number of datapoints whereas 'm' is the number of inputs
[n, m] = size(inputs);
damtx = [];
evs = [];

% 'ind' is an integer which controls the development of new terms 
ind = 1;
greater = 0;
finished = 0;
X = [];
if m == 1
    sett = 1;
elseif way3
    sett = 3;
else
    sett = 2;
end
while 1
    
    indvec = zeros(1,m);
    summ = ind;
            
    while summ
        for j=1:sett
            indvec(j) = indvec(j) + 1;
            summ = summ-1;
            if summ == 0
                break;
            end
        end
    end
    
    while 1
    
        vecs = uniqueperms(indvec);
        
        damtx = [damtx; vecs];
        
        [beters, ~, ~, ~, X, ev] = gibbs(inputs, data, phis, X, damtx, a, b, atau, btau, draws);

        %betavs = [std(beters(ceil(draws/2+1),:))]
        
        if aic
            [dam, ~] = size(damtx);
            ev = ev + (2 - log(n))*dam;
        end
        
        disp([ind ev])
        
        evs = [evs ev];
        if ev == min(evs)
            
            betas = beters;
            mtx = damtx;
            greater = 1;
            
        elseif greater < tolerance
            greater = greater + 1;
        else
            finished = 1;
            break;
        end
        
        if m == 1
            break;
        elseif way3
            if indvec(2) > indvec(3)
                indvec(1) = indvec(1) + 1;
                indvec(2) = indvec(2) - 1;
            elseif indvec(3)
                indvec(2) = indvec(2) + 1;
                indvec(3) = indvec(3) - 1;
                if indvec(2) > indvec(1)
                    indvec(1) = indvec(1) + 1;
                    indvec(2) = indvec(2) - 1;
                end
            else
                break;
            end
        elseif indvec(2)
            indvec(1) = indvec(1) + 1;
            indvec(2) = indvec(2) - 1;
        else
            break;
        end
    
    end
    
    if finished
        break;
    end
    
    ind = ind + 1;
    
    if ind > length(phis)
        break;
    end

end

if gimmie
    betas = beters;
    mtx = damtx;
end

end