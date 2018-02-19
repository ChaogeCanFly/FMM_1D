clc
clear all

N = 15;                %number of Fourier orders
period = 700*(10^(-9));  %period of periodic layer
w = 40*(10^(-9));     %ridge width

%{
L = 2;                 %number of layers
h = zeros(L,1);
%h(1) = 30*10^(-9);       %thickness of periodic layer
%h(2) = 60*10^(-9);
%h(3) = 1*10^(-6);
h(1) = 60*10^(-9);
h(2) = 1*10^(-6);
%}
%gg = gamma(1.7)
%mfun('FresnelC',0:5)
M = 5001;              %number of modes for Fourier transform of epsilon
x = (1:1:M)*period/M;
nnew = 2.028;
epsnew = nnew^2;
refIndices = [1.0 1.0];     

lmin = 740*10^(-9);
lmax = 750*10^(-9);
lambda = linspace(lmin, lmax, 350);
[Nll,Nl] = size(lambda);

theta = [0 0.2 0.4 0.6]*pi/180;
Nt=4;
%{
thetamin = 0;
thetamax=1*pi/180;
theta = linspace(thetamin, thetamax,40);
[Ntt, Nt] = size(theta);
%}

phi = 0*pi/180;
Rsum=zeros(Nl,Nt);
Tsum=zeros(Nl,Nt);

%nTa2O5 = 2.1+0.0001*1j;
%nTa2O5 = 2.1;
%epsTa2O5 = nTa2O5^2;

%nSiO2 = 1.48;
%epsSiO2 = nSiO2^2;
%nSi = 3.71;
%epsilon(:,1)=epsTa2O5*ones(M,1);
%epsilon(:,2)=epsSiO2*ones(M,1);

%lambda = 1.064*(10^(-6));

%theta = 0*pi/180;
%Nt = 1;



%{
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
    Rsum_nonperiodic(i,j) = sum(eta_R1);
    Tsum_nonperiodic(i,j) = sum(eta_T1);
    end
    end
end
%}

L=2;
h(1) = 30*10^(-9);       %thickness of periodic layer
h(2) = 70*10^(-9);
%h(3) = 1*10^(-6);
epsilon = zeros(M, L);
epsilon(:,1)=epsnew*ones(M,1);
epsilon(:,2)=epsnew*ones(M,1);
%epsilon(:,3)=epsSiO2*ones(M,1);
for i=1:M    
    if x(i)<=w     
      epsilon(i,1) = 1.0;
    end
end
eps = zeros(2*N+1,2*N+1,L);
epsinv = zeros(2*N+1,2*N+1,L);
for i=1:L
    eps(:,:,i) = FMM_happy_epsilon_1D(epsilon(:,i), N, M);
    epsinv(:,:,i) = FMM_happy_epsilon_1D(1./epsilon(:,i), N, M);
end
for i=1:Nl
    for j=1:Nt
    for k=1:1
    [beta1, eta_R1, eta_T1] = FMM_1D_TM_RT_multi(eps, epsinv,...
    period, h, lambda(i), theta(j), phi(k), refIndices, N, M, L);
    Rsum(i,j) = sum(eta_R1);
    Tsum(i,j) = sum(eta_T1);
    end
    end
end
%{
for i=1:Nl
 %   Rsum_Fano(i,1)=Rsum(i,1)-Rsum_nonperiodic(i,1);
end

FanoEqn = 'H*(1+2*(x-x0)/(W*q))^2/(1+4*(x-x0)^2/W^2)';
%}
%{
tl = linspace(lmin,lmax,500)*10^9;
tt = linspace(thetamin, thetamax,100)*180/pi;
[XI,YI] = meshgrid(tl,tt);
ZI = griddata(theta*180/pi,lambda*10^9,Rsum,YI,XI);
figure;
%pcolor(kx,ky,Rsum)
pcolor(XI,YI,ZI)
shading flat
caxis([0 1])
%}
%{
Hs = 0.75;
Ws = 0.005;
qs = 1;
x0s = 749.787;
startPoints = [Hs Ws qs x0s];
x = zeros(Nl,1);
y = zeros(Nl,1);
for i=1:Nl
    x(i,1) = lambda(i)*10^9;
    y(i,1) = Rsum_Fano(i);
end
lminnm = lmin*10^9
lmaxnm = lmax*10^9

Rfunction = fit(x, y, FanoEqn, 'start', startPoints, 'Exclude', [lminnm lmaxnm])
Rfit = zeros(Nl,1);
for i=1:Nl
    Rfit(i,1) = Rfunction(x(i));
end 
coeffs = coeffvalues(Rfunction)
figure(6)
hold on
plot(lambda, Rfit, 'r',  'LineWidth', 2)
axis([lmin lmax 0 1])
set(gca,'fontsize', 18)
hold off
%plot(theta*180/pi, Rsum, 'r', 'LineWidth', 2)
%}
%{
figure(1)
hold on
plot(lambda, Rsum(:,1), 'b', 'LineWidth', 2)
%plot(lambda, Rsum_Fano(:,1), 'g', lambda, Rsum(:,1), 'b',  'LineWidth', 2)
h1 = legend('theta=0',1);
set(h1,'Interpreter','none')
axis tight
axis([lmin lmax 0 1])
set(gca,'fontsize', 18)
hold off

figure(2)
hold on
plot(lambda, Rsum(:,2), 'b', 'LineWidth', 2)
h1 = legend('theta=1',1);
set(h1,'Interpreter','none')
axis tight
axis([lmin lmax 0 1])
set(gca,'fontsize', 18)
hold off

figure(3)
hold on
plot(lambda, Rsum(:,3), 'b', 'LineWidth', 2)
h1 = legend('theta=2',1);
set(h1,'Interpreter','none')
axis tight
axis([lmin lmax 0 1])
set(gca,'fontsize', 18)
hold off

figure(4)
hold on
plot(lambda, Rsum(:,4), 'b', 'LineWidth', 2)
h1 = legend('theta=3',1);
set(h1,'Interpreter','none')
axis tight
axis([lmin lmax 0 1])
set(gca,'fontsize', 18)
hold off
%}
figure(5)
hold on
plot(lambda, Rsum(:,1), 'b', lambda, Rsum(:,2)+0.5, 'g', lambda, Rsum(:,3)+1, 'r', lambda, Rsum(:,4)+1.5, 'm', 'LineWidth', 2)
%h5 = legend('theta=0','theta=1','theta=2','theta=3',4);
h5 = legend('theta=0','theta=0.2','theta=0.4','theta=0.6',4);
set(h5,'Interpreter','none')
%title('W=',w*10^9,' nm, D=', period*10^9,' nm')
axis tight
axis([lmin lmax 0 2.5])
set(gca,'fontsize', 18)
hold off

