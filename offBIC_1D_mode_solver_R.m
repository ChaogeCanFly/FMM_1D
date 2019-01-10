clc
clear all

N = 10;                %number of Fourier orders
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

l1 = 800*10^(-9);
dl = -(0.5+0.000j)*10^(-9);
l2 = l1+dl;
lambda = [l1 l2];
c = 3*10^8;
w1 = c/l1;
w2 = c/l2;
dw = w2-w1;
Nl = 2;

kx = 9.06*10^6;
%theta = asin(kx*l1/(2*pi*nmedia));
theta1 = 22*pi/180;
thetamin = theta1;% - 2*pi/180;
thetamax = theta1 + 1*pi/180;
%Nt = 1;
theta = linspace(thetamin,thetamax,10);
[Ntt,Nt] = size(theta);


phi = 0*pi/180;


eps = zeros(2*N+1,2*N+1,L);
epsinv = zeros(2*N+1,2*N+1,L);

for i=1:L
    eps(:,:,i) = FMM_happy_epsilon_1D(epsilon(:,i), N, M);
    epsinv(:,:,i) = FMM_happy_epsilon_1D(1./epsilon(:,i), N, M);
end

N_iterations = 10;

w_eig = zeros(Nt,1);
dw00 = zeros(N_iterations, Nt);
dw0 = w1;

for j=1:Nt
    for ii=1:N_iterations
        for i=1:Nl
            [R] = FMM_1D_TM_multi_mode_solver_R...
                (eps, epsinv, period, h, lambda(i), theta(j), phi, refIndices, N, M, L);
            if i==1
                R0 = R;
            else
                R1 = R;
            end
        end
        dR = (R1-R0)/dw;
        
        A = -dR\R0;
        [V, D] = eig(A);
        
        w_array = sort(diag(D));
        dw0 = min(w_array);
        if ii==1
            dw0 = w_array(1);
        end
        w1 = w1 + dw0;
        w2 = w1 + dw;
        
        l1 = c/w1;
        l2 = c/w2;
        lambda = [l1 l2];
        
        dw00(ii,j) = dw0;
        
    end
    w_eig(j) = w1;
end

lambda_eig = c./w_eig;
k0=2*pi./real(lambda_eig);
%kx = k0*periodx*sin(theta')*cos(phi);
%plot(kx,w_eig,'b','Linewidth',2);
N1 = Nt;
theta_new = theta*180/pi;
%plot(theta_new, lambda_eig*10^9)
plot(lambda_eig*10^9,theta_new)
theta_new1 = theta_new(1:N1);
lambda_eig1 = lambda_eig(1:N1)*10^9;
plot(lambda_eig1*10^9,theta_new1)

p=polyfit(transpose(theta_new1),lambda_eig1,4);
x1 = theta_new1;
y1=polyval(p,x1);
figure(1)
plot(lambda_eig1,theta_new1,'o')
hold on
plot(y1,x1,'r')
hold off

w_eigg=sort(w_eig);
Q=abs(real(w_eigg)./(2*imag(w_eigg)));
figure(2)
plot(theta_new1,Q,'o')
hold off
%{
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
%}
