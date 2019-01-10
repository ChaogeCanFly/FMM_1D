clc
clear all

N = 5;                %number of Fourier orders
period = 1000*(10^(-9));  %period of periodic layer
d = 150*(10^(-9));     %ridge width
h = zeros(3,1);
h(1) = 3*10^(-6);
h(2) = 350*(10^(-9));       %thickness of periodic layer
h(3) = 1*10^(-6);
M = 101;              %number of modes for Fourier transform of epsilon
L = 3;                 %number of layers
x = (1:1:M)*period/M;
epsilon = zeros(M, L);
for i=1:M
    epsilon(i,1) = 1.44^2;
    if x(i)<d
      %  epsilon(i,2) = 3.48^2;
      epsilon(i,2) = 3.48^2;
    else
        epsilon(i,2) = 1.44^2;
    end
    epsilon(i,3) = 1.44^2;
end

%lambda = 1.064*(10^(-6));
refIndices = [1 3.48];     
%{
Rsum=zeros(50,25);
Tsum=zeros(50,25);
lambda = zeros(50,1);
theta = zeros(100,1);
%theta=zeros(90,1);
%theta = 4*pi/180;
%}
lambda = linspace(1400,1500,20)*10^(-9);
[Nll,Nl] = size(lambda)
theta = linspace(1,15,15)*pi/180;
[Ntt, Nt] = size(theta)

theta = [4]*pi/180;
Nt=1;
phi = 0;
Rsum=zeros(Nl,Nt);
Tsum=zeros(Nl,Nt);
for i=1:L
    eps(:,:,i) = FMM_happy_epsilon_1D(epsilon(:,i), N, M);
    epsinv(:,:,i) = FMM_happy_epsilon_1D(1./epsilon(:,i), N, M);
end
for i=1:Nl
    for j=1:Nt
        i
    %lambda(i) = (1200+i*16)*(10^(-9))
    %theta(j) = (j/4)*pi/180;
    [beta1, eta_R1, eta_T1] = FMM_1D_TE_RT_multi(eps, period, h, lambda(i), theta(j), refIndices, N, M, L);
    Rsum(i,j) = sum(eta_R1);
    %Tsum(i) = 1-Rsum(i);
    Tsum(i,j) = sum(eta_T1);
    %eta = cat(2,eta_R1,eta_T1);
    end
end
plot(lambda,Rsum)
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