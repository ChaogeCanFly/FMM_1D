clc
clear all
allPlots = findall(0, 'Type', 'figure', 'FileName', []);
% Close.
delete(allPlots);

N = 5;                  %number of Fourier orders
period = 495*(10^(-9));  %period of periodic layer
w = 150*(10^(-9));        %ridge width
L = 2;                   %number of layers

h = zeros(L,1);
%{
h(1) = 50*10^(-9);       %thickness of periodic layer
h(2) = 40*10^(-9);
h(3) = 1*10^(-6);
%}

h(1) = 90*10^(-9);
h(2) = 1*10^(-6);

M = 5001;              
x = (1:1:M)*period/M;
epsilon = zeros(M, L);
%nTa2O5 = 2.1+0.0001*1j;
nTa2O5 = 2.1;
epsTa2O5 = nTa2O5^2;
nSiO2 = 1.48;
epsSiO2 = nSiO2^2;
nSi = 3.71;
epsilon(:,1)=epsTa2O5*ones(M,1);
%epsilon(:,2)=epsTa2O5*ones(M,1);
epsilon(:,2)=epsSiO2*ones(M,1);
for i=1:M    
    if x(i)<=w     
      epsilon(i,1) = 1.0;
    end
end
refIndices = [1.0 nSi];     

lmin = 600*10^(-9);
lmax = 1000*10^(-9);
lambda = linspace(lmin, lmax, 2500);
[Nll,Nl] = size(lambda);
%theta = 0*pi/180;
%Nt = 1;

%thetamin = 0;
%thetamax=1.0*pi/180;
%theta = linspace(thetamin, thetamax,70);
%[Ntt, Nt] = size(theta);


theta = [0.1 1 3 10]*pi/180;
Nt=4;

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
    end
    end
end

%%%%%%%%%%non-etched%%%%%%%%%%%%%%%%%%%%%%%%%%
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
Rsum = Rsum./Rsum_non;
%{
for i=1:Nl
   Rsum_Fano(i,1)=Rsum(i,1)-Rsum_nonperiodic(i,1);
end

FanoEqn = 'H*(1+2*(x-x0)/(W*q))^2/(1+4*(x-x0)^2/W^2)';
%}
%{
tl = linspace(lmin,lmax,300)*10^9;
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
lambda = lambda*10^6;
lmin = lmin*10^6;
lmax = lmax*10^6;
figure(5)
hold on
plot(lambda, Rsum(:,1), 'b', lambda, Rsum(:,2)+5, 'g', lambda, Rsum(:,3)+10, 'r', lambda, Rsum(:,4)+15, 'm', 'LineWidth', 2)
%plot(lambda, Rsum(:,1), 'b', 'LineWidth', 2)
%h5 = legend('theta=0','theta=1','theta=2','theta=3',4);
h5 = legend('theta=0.1 deg','theta=1 deg','theta=3 deg','theta=10 deg',4);
set(h5,'Interpreter','none')
%title('W=',w*10^9,' nm, D=', period*10^9,' nm')
axis tight
ax = gca;
ax.XAxis.MinorTick = 'on';
%axis([lmin lmax 0 2.5])
%set(gca,'Xtick',0.6:0.025:1)
%grid minor
%set(gca,'XtickLabel',lambda(0.6:0.025:end))
xlabel('lambda, mkm')
ylabel('R')
%set(haxes(1),'XLim',[Min Max],'YLim',[MinV MaxV]);
%Sets Range set(haxes(2),'Xlim',[Min Max],'YLim',[MinI MaxI]);
set(gca,'fontsize', 16)
%xticks([0.6 0.025 1])
%xticks([lmin lmin+0.01 lmax])
%yticks([-1 0 1])
%datetick('x', 'yyyy', 'keepticks');
hold off
