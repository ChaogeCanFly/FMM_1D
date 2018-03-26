clc
clear all
allPlots = findall(0, 'Type', 'figure', 'FileName', []);
% Close.
delete(allPlots);

N = 4;                  %number of Fourier orders
period = 500;  %period of periodic layer
w = 15;       %ridge width
L = 1;         %number of layers

h = zeros(L,1);

h(1) = 100;       %thickness of periodic layer


M = 5001;              
x = (1:1:M)*period/M;
epsilon = zeros(M, L);
%nTa2O5 = 2.1+0.0001*1j;
%nTa2O5 = 2.1;
%epsTa2O5 = nTa2O5^2;
nSi3N4 = 1.8;
epsSi3N4 = nSi3N4^2;
nSiO2 = 1.48;
epsSiO2 = nSiO2^2;
%nSi = 3.71;
epsilon(:,1)=epsSi3N4*ones(M,1);
for i=1:M    
    if x(i)<=w     
      epsilon(i,1) = 1.0;
    end
end
refIndices = [1.0 nSiO2];     

lmin = 740;
lmax = 745;
lambda = linspace(lmin, lmax, 501);
[Nll,Nl] = size(lambda);
%theta = 0*pi/180;
%Nt = 1;

thetamin = 0*pi/180;
thetamax= 0.5*pi/180;
theta = linspace(thetamin, thetamax, 200);
[Ntt, Nt] = size(theta);


%theta = [7 8 9 10 13]*pi/180;
%Nt=5;
%{
theta = [0.1 0.3 1 3]*pi/180;
Nt=4;
%}
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
    [beta1, eta_R1, eta_T1] = FMM_1D_TE_RT_multi(eps, epsinv,...
    period, h, lambda(i), theta(j), phi(k), refIndices, N, M, L);
    Rsum(i,j) = sum(eta_R1);
    Tsum(i,j) = sum(eta_T1);
    lambda(i)
    end
    end
end



%save('sample_497_180_map.mat','theta','lambda','Rsum')

figure(1);
pcolor(lambda,theta*180/pi,transpose(Rsum))
xlabel('lambda, nm');
ylabel('theta, deg');
colormap('jet');
colorbar;
set(gca,'fontsize', 16)
shading flat
%caxis([0 1])
colorbar
hold off
%{
c = physconst('LightSpeed');
h = 4.135666 * 10^(-15);
kx = (2*pi*10^(-6)./lambda)'*sin(theta);
energy = (h*c./lambda');

figure(2);
pcolor(kx,energy,Rsum)
xlabel('kx, mkm^{-1}');
ylabel('energy, eV');
colormap('gray');
colorbar;
set(gca,'fontsize', 16)
shading flat
caxis([0 1])
colorbar
hold off
%}
%{
for i = 1:Nl
    for j=1:Nt
        omega(i) = 2*pi*c/lambda(i);
        k0(i) = 2*pi*1000/lambda(i); %mkm
        kx(i,j) = k0(i)*sin(theta(j));
    end
end
%}
%%%%%%%%%%non-etched%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
Rsum_non=zeros(Nl,Nt);
Tsum_non=zeros(Nl,Nt);
epsilon(:,1)=epsTa2O5*ones(M,1);
eps = zeros(2*N+1,2*N+1,L);
epsinv = zeros(2*N+1,2*N+1,L);
for i=1:L
    eps(:,:,i) = FMM_happy_epsilon_1D(epsilon(:,i), N, M);
    epsinv(:,:,i) = FMM_happy_epsilon_1D(1./epsilon(:,i), N, M);
end
for i=1:Nl
    for j=1:Nt
    for k=1:1
    [beta1, eta_R1, eta_T1] = FMM_1D_TE_RT_multi(eps, epsinv,...
    period, h, lambda(i), theta(j), phi(k), refIndices, N, M, L);
    Rsum_non(i,j) = sum(eta_R1);
    Tsum_non(i,j) = sum(eta_T1);
    end
    end
end
Rsum_normed = Rsum./Rsum_non;
%}


%{
lambda = lambda*10^6;
lmin = lmin*10^6;
lmax = lmax*10^6;
figure(1)
hold on
plot(lambda, Rsum_normed(:,1), 'b', lambda, Rsum_normed(:,2)+5, 'g', lambda, Rsum_normed(:,3)+10,'m',...
    lambda, Rsum_normed(:,4)+15, 'c', lambda, Rsum_normed(:,5)+20, 'k', 'LineWidth', 2)
h5 = legend('theta=7 deg','theta=8 deg','theta=9 deg','theta=10 deg','theta=13 deg', 5);
set(h5,'Interpreter','none')
axis tight
ax = gca;
ax.XAxis.MinorTick = 'on';
xlabel('lambda, mkm')
ylabel('R_normed')
set(gca,'fontsize', 16)
%}
%{
figure(2)
plot(lambda, Rsum(:,1), 'b', lambda, Rsum(:,2)+0.5, 'g', lambda, Rsum(:,3)+1.0,'r',...
    lambda, Rsum(:,4)+1.5, 'm', 'LineWidth', 2)
h5 = legend('theta=0.1 deg','theta=0.3 deg','theta=1 deg','theta=3 deg', 4);
set(h5,'Interpreter','none')
axis tight
ax = gca;
ax.XAxis.MinorTick = 'on';
xlabel('lambda, nm')
ylabel('R')
set(gca,'fontsize', 16)

hold off

figure(3)
plot(lambda,Rsum(:,1),'b','LineWidth', 2)
xlabel('lambda, nm')
ylabel('R')
set(gca,'fontsize', 16)
hold off
%}