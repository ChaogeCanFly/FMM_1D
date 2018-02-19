function [beta, W, pplus, pminus] = FMM_1D_TM_RT_beta_e(eps, epsinv, period, h, lambda, theta, refIndices, N, M)
%
% function for one layer of periodic structure
% 
% OUTPUT:
% <beta>: vector with propagation constants of grating modes
% <e>: electric field vectors
% <W>: interface matrix
% <Pl, Pr>: left and right layer matrix
% program for TE polarization
% for TM polarization we switch eps <-> -mu and E<->H

k0 = 2*pi/lambda;
kxv = zeros(2*N+1,1);
for m=1:(2*N+1)
    kxv(m) = (refIndices(1)*k0*sin(theta)+(m-N-1)*2*pi/period);
end
kx = diag(kxv);

miden = eye(2*N+1,2*N+1);
matr = epsinv\((k0^2)*miden - kx/eps*kx);
[hh, beta2] = eig(matr);
beta = sqrt(beta2);
%ff = real(beta)+imag(beta)

pplusv = zeros(2*N+1,1);
pminusv = zeros(2*N+1,1);

for m=1:(2*N+1)
    pplusv(m) = exp(1i*beta(m,m)*h);
    pminusv(m) = exp(-1i*beta(m,m)*h);
end;
pplus = diag(pplusv);
pminus = diag(pminusv);
mzero = zeros(2*N+1, 2*N+1);
miden = eye(2*N+1, 2*N+1);

Pl_up = cat(2, miden, mzero);
Pl_down = cat(2, mzero, pplus);
Pl = cat(1, Pl_up, Pl_down);

Pr_up = cat(2,pplus, mzero);
Pr_down = cat(2, mzero, miden);
Pr = cat(1, Pr_up, Pr_down);

%p_up = cat(2, pplus, mzero);
%p_down = cat(2, mzero, pminus);
%P = cat(1, p_up, p_down);

e = 1/k0*(epsinv*hh*beta);
W_up = cat(2, hh, hh);
W_down = cat(2, e, -e);
W = cat(1, W_up, W_down);

%eigenvalues beta - kz in periodic layer
%eigenvectors e - electric field amplitudes E
%each column of e consists of Fourier coefficients with m from -N to N
%each column corresponds to a Rayleigh mode with different q
end