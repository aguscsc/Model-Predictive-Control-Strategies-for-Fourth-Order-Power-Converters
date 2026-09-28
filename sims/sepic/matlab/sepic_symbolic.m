%% Symbolic characteristic polynomial of the SEPIC small-signal model
syms s D oD L1 L2 Cs Co R positive
% D = duty cycle, oD = D' = 1 - D (kept as a separate symbol for readability)

A_sym = [0            0            -oD/L1   -oD/L1;
         0            0             D/L2    -oD/L2;
         oD/Cs  -D/Cs     0            0;
         oD/Co    oD/Co     0           -1/(R*Co)];

% Characteristic polynomial, grouped by powers of s
charpoly_s = collect(expand(det(s*eye(4) - A_sym)), s)

% Individual coefficients, highest power first: [1; a3; a2; a1; a0]
a_sym = coeffs(charpoly_s, s, 'All').';
a_sym = simplify(subs(a_sym, D, 1 - oD))   % applies D + D' = 1
pretty(a_sym)
pretty(charpoly_s)


%% Factorization check at rho = 1  (L1 = D'/k, L2 = D/k)
syms k positive
cp_rho1 = subs(charpoly_s, [L1, L2], [oD/k, D/k]);
cp_rho1 = simplify(subs(cp_rho1, D, 1 - oD));
factor(cp_rho1, s)