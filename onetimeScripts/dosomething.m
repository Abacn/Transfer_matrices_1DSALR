% print secondB and thirdB to a file

coeffs = [1, 2.5, 4, 1];
xis=[0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1, 2, 3, 4, 5];
Ts = linspace(0.06, 0.4, 18);
for xi=xis
    coeff = [coeffs, xi];
    sb = secondB(
end