function [Rtotal] = FMM_1D_TM_multi_mode_solver_R(eps, epsinv,...
    period, h, lambda, theta, phi, refIndices, N, M, L)
%
% INPUT:
% <epsilon>: vector with permitivity distribution
% <period>: period of the grating in micron
% <lambda>: wavelength in micron
% <theta>: angle of incidence in degree
% <refIndices>: vector [n1 n2] with refr.indices of surrounding
% <N>: number of Fourier orders


% OUTPUT:
% <beta>: vector with propagation constants of grating modes
% program for TE polarization
% for TM polarization we switch eps <-> -mu and E<->H

%after Li,2003 article

k0 = 2*pi/lambda;
beta = zeros(2*N+1,2*N+1,L);

W = zeros(4*N+2,4*N+2,L);
pplus = zeros(2*N+1,2*N+1,L);
pminus = zeros(2*N+1,2*N+1,L);

for i=1:L
    [betatemp,Wtemp,pplust,pminust] = FMM_1D_TM_RT_beta_e(eps(:,:,i), epsinv(:,:,i), period, h(i), lambda, theta, refIndices, N, M);
    beta(:,:,i) = betatemp;
    W(:,:,i) = Wtemp;
    pplus(:,:,i) = pplust;
    pminus(:,:,i) = pminust;
end


kxv = zeros(2*N+1,1);
for m=1:(2*N+1)
    kxv(m) = (refIndices(1)*k0*sin(theta)+(m-N-1)*2*pi/period);
end
kx = diag(kxv);
kz1v = zeros(2*N+1,1);
kz2v = zeros(2*N+1,1);

for m=1:(2*N+1)
    kz1v(m) = ( (refIndices(1)^2)*(k0^2) - (kxv(m))^2 )^(1/2);
    kz2v(m) = ( (refIndices(2)^2)*(k0^2) - (kxv(m))^2 )^(1/2);
end;
mzero = zeros(2*N+1, 2*N+1);
miden = eye(2*N+1, 2*N+1);
kz1 = diag(kz1v);
kz2 = diag(kz2v);
kz1_TM = kz1/(k0*(refIndices(1)^2));
kz2_TM = kz2/(k0*(refIndices(2)^2));

K_up = cat(2, miden, miden);
K1_down = cat(2, kz1_TM, - kz1_TM);
K1 = cat(1, K1_down, K_up);
K2_down = cat(2, kz2_TM, - kz2_TM);
K2 = cat(1, K2_down, K_up);

Smin1 = eye(4*N+2,4*N+2);
S0 = new_recursion(Smin1, K1, W(:,:,1), eye(2*N+1,2*N+1),eye(2*N+1,2*N+1),N);
Stemp = S0;

for i=1:(L-1)
    Si = new_recursion(Stemp, W(:,:,i), W(:,:,i+1), pplus(:,:,i), pminus(:,:,i), N);
    Stemp = Si;
end

Rtotal = new_recursion(Stemp, W(:,:,L), K2, pplus(:,:,L), pminus(:,:,L), N);
