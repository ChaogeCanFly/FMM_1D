function [beta, eta_R, eta_T] = FMM_1D_TE_RT_multi(eps, epsinv, period, h, lambda, theta, phi, refIndices, N, M, L)
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
    [betatemp,Wtemp,pplust,pminust] = FMM_1D_TE_RT_beta_e(eps(:,:,i), period, h(i), lambda, theta, refIndices, N, M);
    beta(:,:,i) = betatemp;
    W(:,:,i) = Wtemp;
    pplus(:,:,i) = pplust;
    pminus(:,:,i) = pminust;
end


kxv = zeros(2*N+1,1);
for m=1:(2*N+1)
    kxv(m) = (refIndices(1)*k0*sin(theta)*cos(phi)+(m-N-1)*2*pi/period);
end
ky = refIndices(1)*k0*sin(theta)*sin(phi);
kx = diag(kxv);
kz1v = zeros(2*N+1,1);
kz2v = zeros(2*N+1,1);

for m=1:(2*N+1)
    kz1v(m) = ( (refIndices(1)^2)*(k0^2) - (kxv(m))^2 - ky^2 )^(1/2);
    kz2v(m) = ( (refIndices(2)^2)*(k0^2) - (kxv(m))^2 - ky^2 )^(1/2);
end;
mzero = zeros(2*N+1, 2*N+1);
miden = eye(2*N+1, 2*N+1);
kz1 = diag(kz1v);
kz2 = diag(kz2v);

K_up = cat(2, miden, miden);
K1_down = cat(2, kz1, - kz1);
K1 = cat(1, K_up, K1_down);
K2_down = cat(2, kz2, - kz2);
K2 = cat(1, K_up, K2_down);

t=zeros(4*N+2,4*N+2,L-1);
s_int=zeros(4*N+2,4*N+2,L-1);
s_layer=zeros(4*N+2,4*N+2,L-1);


Smin1 = eye(4*N+2,4*N+2);
S0 = new_recursion(Smin1, K1, W(:,:,1), eye(2*N+1,2*N+1),eye(2*N+1,2*N+1),N);
Stemp = S0;

for i=1:(L-1)
    Si = new_recursion(Stemp, W(:,:,i), W(:,:,i+1), pplus(:,:,i), pminus(:,:,i), N);
    Stemp = Si;
end

Stotal = new_recursion(Stemp, W(:,:,L), K2, pplus(:,:,L), pminus(:,:,L), N);

%{
Smin1 = eye(4*N+2,4*N+2);
S0 = new_recursion_new(Smin1, K2, W(:,:,1), eye(2*N+1,2*N+1), pminus(:,:,1), N);
Stemp = S0;
if L>1
    for i=1:(L-1)
        Si = new_recursion_new(Stemp, W(:,:,i), W(:,:,i+1), pplus(:,:,i), pminus(:,:,i+1), N);
        Stemp = Si;
    end
end
Stotal = new_recursion_new(Stemp, W(:,:,L), K1, pplus(:,:,L), eye(2*N+1,2*N+1), N);
%}
u0 = zeros(2*N+1,1);
u0(N+1) = 1;
%{
u0(N+1) = cos(phi);
Ey_0= u0(N+1);
Ex_0= -Ey_0*tan(phi);
Ez_0= -(kxv(N+1)*Ex_0+ky*Ey_0)/kz1v(N+1);
norm  = Ex_0^2+Ey_0^2+Ez_0^2;
u0=u0/sqrt(norm);
%}
%norm = u0(N+1)^2;
%u0=u0/norm;
dlast = zeros(2*N+1,1);
%dlast(N+1) = 1;
del = cat(1,u0,dlast);

T_R = Stotal*del;
T = T_R(1:(2*N+1));        %u_last
R = T_R((2*N+2):(4*N+2));  %d0
eta_R = zeros(2*N+1,1);
eta_T = zeros(2*N+1,1);
Ex_r=zeros(2*N+1,1);
Ey_r=zeros(2*N+1,1);
Ez_r=zeros(2*N+1,1);
Ex_t=zeros(2*N+1,1);
Ey_t=zeros(2*N+1,1);
Ez_t=zeros(2*N+1,1);



R_new=zeros(2*N+1,1);
T_new=zeros(2*N+1,1);

for m=1:(2*N+1)
%{    
    phi(m) = atan(ky/kxv(m))
    Ey_r(m)= R(m);
    Ex_r(m)= -Ey_r(m)*tan(phi(m));
    Ez_r(m)= -(kxv(m)*Ex_r(m)+ky*Ey_r(m))/kz1v(m);
    
    Ey_t(m)= T(m);
    Ex_t(m)= -Ey_t(m)*tan(phi(m));
    Ez_t(m)= -(kxv(m)*Ex_t(m)+ky*Ey_t(m))/kz2v(m);
    
    R_new(m) = sqrt(Ex_r(m)^2+Ey_r(m)^2+Ez_r(m)^2);
    T_new(m) = sqrt(Ex_t(m)^2+Ey_t(m)^2+Ez_t(m)^2);
    eta_R(m) = (abs(R_new(m)))^2*real( kz1v(m)/(refIndices(1)*k0*cos(theta)) );
    eta_T(m) = (abs(T_new(m)))^2*real( kz2v(m)/(refIndices(2)*k0*cos(theta)) );
%}    
    eta_R(m) = (abs(R(m))^2)*real( kz1v(m)/(refIndices(1)*k0*cos(theta)) );
    eta_T(m) = (abs(T(m))^2)*real( kz2v(m)/(refIndices(2)*k0*cos(theta)) );
   
end
%ex0 = Ex_r(N+1)
%ey0 = Ey_r(N+1)
%ez0 = Ez_r(N+1)

end