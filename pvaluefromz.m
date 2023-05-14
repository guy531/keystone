%% Return the exact p value (one tail and 2 tails), from a given z (standardized distribution)
% INPUT:
%       z value 
% OUTPUT:
%       one tail and two tails exact p value for that z
function [pvalue] =  pvaluefromz (z)


    z = abs(z);
    F = @(x)(exp (-0.5*(x.^2))./sqrt (2*pi));
    p = quad (F, z, 100);
    fprintf ('\nOne tail p-value : %1.6f', p);

    fprintf ('\nTwo tail p-value : %1.6f\n' ,p*2)

    %Return two tail p-value
%     pvalue = p*2; 
    pvalue = p;
    
end