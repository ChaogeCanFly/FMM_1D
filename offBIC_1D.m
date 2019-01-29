clc
clear all

N = 10;                %number of Fourier orders
L = 1;                 %number of layers
h = zeros(L,1);

period = 382*10^(-9);  %period of periodic layer
th = 2.5*period;
h(1) = th;
w = 0.45*period;     %ridge width

M = 5001;              %number of modes for Fourier transform of epsilon
x = (1:1:M)*period/M;
epsilon = zeros(M, L);
nlattice = 2.125;
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

wmin = 350*10^12;
wmax = 270*10^12;
c = 3*10^8;
lmin = c/wmin;
lmax = c/wmax;
lambda = linspace(lmin, lmax, 50);
[Nll,Nl] = size(lambda);

thetamin = 19*pi/180;
thetamax = 89*pi/180;
theta = linspace(thetamin, thetamax,50);
[Ntt, Nt] = size(theta);
%theta=12*pi/180;
%Nt=1;

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
%{
figure;
plot(lambda,Rsum,'r','Linewidth',2)
hold off
%}
for i=1:Nl
    for j=1:Nt
        if Rsum(i,j)>
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
