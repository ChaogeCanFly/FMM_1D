function [beta, eta_R, eta_T] = FMM_1D_TE_RT_book(epsilon, period, h, lambda, theta, refIndices, N, M)
%
% INPUT:
% <epsilon>: vector with permitivity distribution
% <period>: period of the grating in micron
% <lambda>: wavelength in micron
% <theta>: angle of incidence in degree
% <refIndices>: vector [n1 n2] with refr.indices of surrounding
% <N>: number of Fourier orders
%
% OUTPUT:
% <beta>: vector with propagation constants of grating modes
% program for TE polarization
% for TM polarization we switch eps <-> -mu and E<->H

k0 = 2*pi/lambda;

eps_m = fft(epsilon)/M;
e1 = eps_m(1:(2*N+1));
e2 = eps_m((M-2*N):M);

%we have to renumerate e2 for Toeplitz matrix, because matlab
%writes eps(m) with negative m in reverse order
for i=(1:N)
    b = e2(i+1);
    e2(i+1) = e2(2*N+1-i+1);
    e2(2*N+1-i+1)=b;
end;
e2(1)=e1(1);

eps = toeplitz(e1,e2);
kxv = zeros(2*N+1,1);
for m=1:(2*N+1);
    kxv(m) = (refIndices(1)*k0*sin(theta)+(m-N-1)*2*pi/period);
end
kx = diag(kxv);
%n1 = size(eps)
%n2 = size(kx2)
matr = (k0^2)*eps - kx^2;
%eigenvalues beta - kz in periodic layer
%eigenvectors e - electric field amplitudes E
%each column of e consists of Fourier coefficients with m from -N to N
%each column corresponds to a Rayleigh mode with different q
[e, beta2] = eig(matr);
beta = sqrt(beta2);


kz1v = zeros(2*N+1,1);
kz2v = zeros(2*N+1,1);
pplusv = zeros(2*N+1,1);
pminusv = zeros(2*N+1,1);

for m=1:(2*N+1)
    kz1v(m) = ( (refIndices(1)^2)*(k0^2) - (kxv(m))^2 )^(1/2);
    kz2v(m) = ( (refIndices(2)^2)*(k0^2) - (kxv(m))^2 )^(1/2);
    pplusv(m) = exp(1i*beta(m,m)*h);
    pminusv(m) = exp(-1i*beta(m,m)*h);
end;
mzero = zeros(2*N+1, 2*N+1);
miden = eye(2*N+1, 2*N+1);
kz1 = diag(kz1v);
kz2 = diag(kz2v);
pplus = diag(pplusv);
pminus = diag(pminusv);


K_up = cat(2, miden, miden);
K1_down = cat(2, kz1, - kz1);
K1 = cat(1, K_up, K1_down);
K2_down = cat(2, kz2, - kz2);
K2 = cat(1, K_up, K2_down);

e_beta = e*beta;
E_up = cat(2, e, e);
E_down = cat(2, e_beta, - e_beta);
E = cat(1, E_up, E_down);

p_up = cat(2, pplus, mzero);
p_down = cat(2, mzero, pminus);
P = cat(1, p_up, p_down);

t = (K2\E)*(P/E)*K1;
%t=K2\K1; %for two medias only
t11 = t( (1:(2*N+1)), (1:(2*N+1)));
t21 = t( ((2*N+2):(4*N+2)), (1:(2*N+1)) );
t12 = t( (1:(2*N+1)), ((2*N+2):(4*N+2)) );
t22 = t( ((2*N+2):(4*N+2)), ((2*N+2):(4*N+2)) );

%S-matrix
s11 = t11-(t12/t22)*t21;
s12 = t12/t22;
s21 = -t22\t21;
s22 = miden/t22;
s_up = cat(2,s11,s12);
s_down = cat(2,s21,s22);
s=cat(1,s_up,s_down);

deltam0 = zeros(2*N+1,1);
deltam0(N+1) = 1;
delta0 = zeros(2*N+1,1);
del = cat(1,deltam0,delta0);

new_t_up = cat(2, miden, -t12);
new_t_down = cat(2, mzero, -t22);
new_t = cat(1, new_t_up, new_t_down);
t_right = cat(1, t11, t21);
T_R = (new_t\t_right)*deltam0;
%T_R = s*del;
T = T_R(1:(2*N+1));
R = T_R((2*N+2):(4*N+2));
eta_R = zeros(2*N+1,1);
eta_T = zeros(2*N+1,1);
for m=1:(2*N+1)
    %eta_R(m) = (abs(R(m))^2);
    %strictly speaking, in Fresnel's formulas eta_T=1-eta_R,
    %eta_T is not defined as the ratio of the transmitted wave’s
    %electric field amplitude to that of the incident wave,
    %because waves with the same amplitude have different energy
    %in different medias
    
    %diffraction efficiencies are
    eta_R(m) = (abs(R(m))^2)*real( kz1v(m)/(refIndices(1)*k0*cos(theta)) );
    eta_T(m) = (abs(T(m))^2)*real( kz2v(m)/(refIndices(2)*k0*cos(theta)) );
    %The diffraction efficiency of a propagating order is defined
    %as the ratio of the z component of
    %the Poynting vector of that order to the z component
    %of the Poynting vector of the incident plane wave.
end

end