clc
clear all

N = 5;                %number of Fourier orders
period = 500*(10^(-9));  %period of periodic layer
w = 30*(10^(-9));     %ridge width
h = zeros(1,1);
h(1) = 90*10^(-9);
%h(2) = 90*(10^(-9));       %thickness of periodic layer
%h(1) = 1*10^(-6);
M = 101;              %number of modes for Fourier transform of epsilon
L = 1;                 %number of layers
x = (1:1:M)*period/M;
epsilon = zeros(M, L);
for i=1:M
    %epsilon(i,1) = 1.48^2;
    if x(i)<w     
      epsilon(i,1) = 1.0;
    else
        epsilon(i,1) = 2.1^2;
    end
end

%lambda = 1.064*(10^(-6));
refIndices = [1.0 1.48];     

%lambda = linspace(1400,1500,50)*10^(-9);
%[Nll,Nl] = size(lambda)
lambda = 780*10^(-9);
Nl=1;
theta = linspace(0,89,89)*pi/180;
[Ntt, Nt] = size(theta)


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
    
    [beta1, eta_R1, eta_T1] = FMM_1D_TE_RT_multi(eps, epsinv, period, h, lambda(i), theta(j), phi(k), refIndices, N, M, L);
    Rsum(i,j) = sum(eta_R1);
    %Tsum(i) = 1-Rsum(i);
    Tsum(i,j) = sum(eta_T1);
    %eta = cat(2,eta_R1,eta_T1);
        end
    end
end

figure(2)
hold on
plot(theta*180/pi, Rsum, 'r', 'LineWidth', 2)
set(gca,'fontsize', 18)
hold off
%{
figure;
imagesc(theta,lambda,Rsum);
set(gca,'Yscale','linear','Ydir','normal');
colormap gray;
%}
%num = zeros(2*N+1,1);
%for q=1:(2*N+1)
%   num(q) = q-N-1;
%end
%bar(num,eta,'stack')
%bar(num, eta_T1)