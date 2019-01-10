clc
clear all

N = 8;                %number of Fourier orders
L = 1;                 %number of layers
h = zeros(L,1);
th = 212*10^(-9);
h(1) = th;
period = 1.8*th;  %period of periodic layer
w = 0.45*period;     %ridge width

M = 5001;              %number of modes for Fourier transform of epsilon
x = (1:1:M)*period/M;
epsilon = zeros(M, L);
nlattice = 2.08;
epslattice = nlattice^2;
nmedia = 1.46;
epsmedia = nmedia^2;

refIndices = [nmedia nmedia];
epsilon(:,1)=epsmedia*ones(M,1);

for i=1:M
    if x(i)<w
        epsilon(i,1) = epslattice;
    end
end

lmin = 710*10^(-9);
lmax = 724*10^(-9);
lambda = linspace(lmin, lmax, 56);
[Nll,Nl] = size(lambda);
%theta = 0*pi/180;
%Nt = 1;
thetamin = 12*pi/180;
thetamax = 15*pi/180;
theta = linspace(thetamin, thetamax,60);
[Ntt, Nt] = size(theta);


phi = 0*pi/180;
Rsum=zeros(Nl,Nt);
Tsum=zeros(Nl,Nt);

eps = zeros(2*N+1,2*N+1,L);
epsinv = zeros(2*N+1,2*N+1,L);

for i=1:L
    eps(:,:,i) = FMM_happy_epsilon_1D(epsilon(:,i), N, M);
    epsinv(:,:,i) = FMM_happy_epsilon_1D(1./epsilon(:,i), N, M);
end
for i=1:Nl
    for j=1:Nt
    for k=1:1
    [beta1, eta_R1, eta_T1] = FMM_1D_TM_RT_multi(eps, epsinv, period, h, lambda(i), theta(j), phi(k), refIndices, N, M, L);
    Rsum(i,j) = sum(eta_R1);
    Tsum(i,j) = sum(eta_T1);
    end
    end
end
XI = lambda;
YI = theta*180/pi;
ZI = transpose(Rsum);
figure;
pcolor(XI,YI,ZI)
shading flat
colormap('jet');
colorbar;
%caxis([0 0.4])
hold off

