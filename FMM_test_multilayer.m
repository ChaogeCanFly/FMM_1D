clc
clear all

N = 15;                %number of Fourier orders
period = 5*(10^(-6));  %period of periodic layer
d = 2.5*(10^(-6));     %ridge width
h = 1*(10^(-6));       %thickness of periodic layer
M = 1001;              %number of modes for Fourier transform of epsilon
L = 2;                 %number of layers
x = (1:1:M)*period/M;
epsilon = zeros(M, L);
for i=1:M
    if x(i)<d
        epsilon(i,1) = 9.0;
    else
        epsilon(i,1) = 9.0;
    end
    epsilon(i,2) = 9.0;
end

lambda = 1.064*(10^(-6));
refIndices = [3.0 1.0];     

Rsum=zeros(90,1);
Tsum=zeros(90,1);
theta=zeros(90,1);

for i=1:1:89
    theta(i) = (i-1)*pi/180;
    [beta1, eta_R1, eta_T1] = FMM_1D_TE_RT_multi(epsilon, period, h, lambda, theta(i), refIndices, N, M, L);
    Rsum(i) = sum(eta_R1);
    %Tsum(i) = 1-Rsum(i);
    Tsum(i) = sum(eta_T1);
    %eta = cat(2,eta_R1,eta_T1);
end

plot([0:1:89],Rsum,'b',[0:1:89],Tsum,'r')

%num = zeros(2*N+1,1);
%for q=1:(2*N+1)
%   num(q) = q-N-1;
%end
%bar(num,eta,'stack')
%bar(num, eta_T1)