
function [meen, bounds, rmse] = coverage(betas, normputs, data, phis, mtx, draws, plots)

[m, mbets] = size(betas);

[n,mputs] = size(normputs);

setnos_p = randi(m, [1 draws]);
while (1)
    setnos = unique(setnos_p);
    if length(setnos) == length(setnos_p)
        break;
    else
        setnos_p = [setnos randi(m,[1 draws-length(setnos)])];
    end
end

X = zeros(n, mbets);

for i=1:n

    phind = ceil(normputs(i,:)*499);
    phind = phind + (phind == 0);

    for j=1:mbets-1

        phi = 1;
        for k=1:mputs

            num = mtx(j,k);

            if num
                x = 499*normputs(i,k) - phind(k) + 1;
                phi = phi*(phis{num}.zero(phind(k)) + phis{num}.one(phind(k))*x + phis{num}.two(phind(k))*x^2 + phis{num}.three(phind(k))*x^3);
            end
        end

        X(i,j) = phi;

    end
end

X(:,mbets) = ones(n,1);

modells = zeros(length(data), draws);
for i=1:draws
    modells(:,i) = X*betas(setnos(i), :)';
end

meen = mean(modells,2);
bounds = zeros(length(data),2);
cut = floor(draws*0.025);
for i=1:length(data)
    drawset = sort(modells(i,:));
    bounds(i,1) = drawset(cut);
    bounds(i,2) = drawset(draws-cut);
end


if (plots)
    figure
    hold on

    plot(meen, 'b', 'LineWidth', 2.0)
    plot(bounds(:,1), 'k--')
    plot(bounds(:,2), 'k--')

    plot(data, 'ro')
end

rmse = sqrt(mean((meen-data).^2));

end