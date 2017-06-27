% SALR boltzmanweight
function wei = boltzmanweight(rs, coeffs)
    % coeffs = [\sigma, \lambda, \kappa, \beta*\epsilon, \xi]
    sigma = coeffs(1);
    wei = zeros(size(rs));
    index = 1;
    for r=abs(rs(:)).'
        if(r<sigma) 
            u = 0;
        elseif (r<coeffs(2)*sigma)
            u = exp(coeffs(4));
        elseif (r<coeffs(3)*sigma)
            u = exp(-coeffs(5)*coeffs(4)*(coeffs(3)-r/sigma));
        else
            u = 1;
        end
        wei(index) = u;
        index = index + 1;
    end;
end
