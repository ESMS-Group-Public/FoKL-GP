
function phi = splineconvert500(A)

% for Yinkai's orthogonal basis functions that are based on a normalized
% interval

% dx = 1/499;
% x = 0:dx:1-dx;
% x = x';

%wunn = ones(size(x));

phi = cell(500,1);

for i=1:500
    
    a = A((i-1)*499+1:i*499,4);
    b = A((i-1)*499+1:i*499,3);
    c = A((i-1)*499+1:i*499,2);
    d = A((i-1)*499+1:i*499,1);

    phi{i}.three = a; %a/dx^3;
    phi{i}.two = b; %b.*(wunn - 3*a.*x.^3/dx^3)/dx^2;
    phi{i}.one = c; %c.*(wunn + 3*a.*x.^2/dx^2 - 2*b.*x/dx^2)/dx;
    phi{i}.zero = d; %d - a.*x.^3/dx^3 + b.*x.^2/dx^2 - c.*x/dx;
    
end

end