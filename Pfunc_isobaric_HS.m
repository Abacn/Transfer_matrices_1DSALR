% Calculate the partition function, 1D hard sphere
function pfunc = Pfunc_isobaric_HS(P, beta, coeffs, divide)
    % set 2*pi*m/h^2 = 1
    pfunc = exp(-beta*P*coeffs(1))./(beta*P);
    pfunc = sqrt(1/beta)*pfunc;
end