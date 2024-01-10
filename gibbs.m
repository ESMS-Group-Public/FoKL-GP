classdef gibbs
    %sampler class

    properties
        phis 
        inputs
        ranges
        data
        a
        b
        atau
        btau
        chimod
        mh
        basisfile
        modparams
    end

    methods
        function obj = gibbs(inputsin,datain,ain,bin,atauin,btauin,chimodin,modparamsin,rangesin,basisin)

            obj.inputs = inputsin;
            obj.data = datain;
            obj.a = ain;
            obj.b = bin;
            obj.atau = atauin;
            obj.btau = btauin;
            obj.chimod = chimodin;
            obj.ranges = rangesin;
            obj.modparams = modparamsin;

            obj.basisfile = basisin;
            obj = splineconvert500(obj);
        end

        function [betas, sigs, taus, betahat, X, ev, obj] = sampler(obj, Xin, discmtx, draws)

            % this version of the sampler increases efficiency by accepting an set of
            % inputs 'Xin' derived from earlier iterations.

            % 'inputs' is the set of normalized inputs -- both parameters and model inputs -- with
            % columns corresponding to inputs and rows the different
            % experimental designs

            % 'data' are the experimental results: column vector, with entries
            % corresponding to rows of 'inputs'

            % 'chimod' is a string specifying the function handle to be used to evaluate X

            % 'discmtx' is the interaction matrix for the bss-anova function -- rows
            % are terms in the function and columns are inputs (cols should line up
            % with cols in 'inputs'

            % 'a' and 'b' are the parameters of the ig distribution for the
            % observation error variance of the data

            % 'atau' and 'btau' are the parameters of the ig distribution for the 'tau
            % squared' parameter: the variance of the beta priors

            % 'draws' is the total number of draws


            %% first build the matrix by calculating the corresponding basis function outputs for each set of inputs
            [minp, ninp] = size(obj.inputs);

            [mmtx, ~] = size(discmtx);
            if isempty(Xin)
                Xin = ones(minp,1);
            end
            [~, nxin] = size(Xin);

            if isequal(obj.chimod,'standard')
                X = [Xin zeros(minp, mmtx-nxin+1)]; % number of data points by number of terms in the function

                for i=1:minp

                    phind = ceil(obj.inputs(i,:)*499);
                    phind = phind + (phind == 0);

                    for j=nxin:mmtx

                        phi = 1;
                        for k=1:ninp

                            num = discmtx(j,k);

                            if num
                                x = 499*obj.inputs(i,k) - phind(k) + 1;
                                phi = phi*(obj.phis{num}.zero(phind(k)) + obj.phis{num}.one(phind(k))*x + obj.phis{num}.two(phind(k))*x^2 + obj.phis{num}.three(phind(k))*x^3);
                            end
                        end

                        X(i,j+1) = phi;

                    end
                end
            elseif isequal(obj.chimod, 'standardC')
                % persistent mh
                % if ~(isa(mh,'matlab.mex.MexHost') && isvalid(mh))
                %     mh = mexhost;
                % end
                obj.mh = mexhost;
                Xapp = feval(obj.mh, 'chimatrix_eval', obj.inputs, discmtx(nxin:mmtx,:));
                X = [Xin Xapp(:,2:end)];
            else
                Xapp = feval(obj.chimod, obj.inputs, obj.ranges, obj.modparams, obj.phis, discmtx(nxin:mmtx,:));
                Xin(:,1) = Xapp(:,1);
                X = [Xin Xapp(:,2:end)];
            end

            % initialize sigsqd and tausqd at the mode of their priors
            sigsqd = obj.b/(1+obj.a);
            tausqd = obj.btau/(1+obj.atau);

            XtX = X'*X;

            Xty = X'*obj.data;

            [Q, Lamb] = eig(XtX);
            QtXty = Q'*Xty;

            Lamb_inv = diag(1./diag(Lamb));

            betahat = Q*Lamb_inv*QtXty;
            %squerr = norm(data - X*betahat)^2;

            astar = obj.a + length(obj.data)/2 + (mmtx + 1)/2;
            atau_star = obj.atau + (mmtx+1)/2;

            dtd = obj.data'*obj.data;

            %% Gibbs iterations

            betas = zeros(draws, mmtx+1);
            sigs = zeros(draws, 1);
            taus = zeros(draws, 1);

            %lik = zeros(draws, 1);
            n = length(obj.data);

            for k = 1:draws

                Lamb_tausqd = Lamb + (1/tausqd)*eye(mmtx+1);
                Lamb_tausqd_inv = diag(1./diag(Lamb_tausqd));

                mun = Q*Lamb_tausqd_inv*QtXty;
                S = Q*diag(sqrt(diag(Lamb_tausqd_inv)));

                vec = normrnd(0,1,[mmtx+1,1]);
                betas(k,:) = (mun + sqrt(sigsqd)*S*vec)';

                %lik(k) = -(n/2)*log(sigsqd) - (squerr + (betahat' - betas(k,:))*XtX*(betahat - betas(k,:)'))/(2*sigsqd);
                %lik(k) = -(n/2)*log(sigsqd) - norm(data - X*betas(k,:)')^2/(2*sigsqd);
                %evs(k) = sqrt(prod(diag(Lamb_tausqd)))/(2*pi*sigsqd)^(mmtx/2 + 0.5);

                bstar = obj.b + 0.5*(betas(k,:)*XtX*betas(k,:)' - 2*betas(k,:)*Xty + dtd + betas(k,:)*betas(k,:)'/tausqd);
                sigsqd = 1/gamrnd(astar, 1/bstar);

                sigs(k) = sigsqd;
                btau_star = (1/(2*sigsqd))*(betas(k,:)*betas(k,:)') + obj.btau;
                tausqd = 1/gamrnd(atau_star, 1/btau_star);
                taus(k) = tausqd;

            end

            %% calculate the evidence

            siglik = var(obj.data - X*betahat);
            lik = -(n/2)*log(siglik) - (n-1)/2;
            ev = (mmtx+1)*log(n) - 2*max(lik);

        end

    end
    methods (Access = protected)    
        function obj = splineconvert500(obj)

            % spline coefficients based on a normalized
            % interval

            A = readmatrix(obj.basisfile);
            obj.phis = cell(500,1);

            for i=1:500

                aa = A((i-1)*499+1:i*499,4);
                bb = A((i-1)*499+1:i*499,3);
                c = A((i-1)*499+1:i*499,2);
                d = A((i-1)*499+1:i*499,1);

                obj.phis{i}.three = aa; %a/dx^3;
                obj.phis{i}.two = bb; %b.*(wunn - 3*a.*x.^3/dx^3)/dx^2;
                obj.phis{i}.one = c; %c.*(wunn + 3*a.*x.^2/dx^2 - 2*b.*x/dx^2)/dx;
                obj.phis{i}.zero = d; %d - a.*x.^3/dx^3 + b.*x.^2/dx^2 - c.*x/dx;

            end

        end

    end
end