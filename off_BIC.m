
clc;
clear all;
%% initializations
Nt = 50;
Nl = 50;
wmin = 350*10^12;
wmax = 270*10^12;
c = 3*10^8;
lmin = c/wmin*10^6;
lmax = c/wmax*10^6;
wll = linspace(0.71,0.724,Nl); % wavelength
theta_min = 19*pi/180;
theta_max = 89*pi/180;
theta = linspace(theta_min,theta_max,Nt);
kxx0 = zeros(Nt,1);
for j=1:Nt
    kxx0(j) = sin(theta(j));
end
% incidence wavevector horizontal projection (dimensionless) (Bloch wavevector)

gh = 0.212; % grating depth
gp = 1.8*gh; % grating period
no = 100; % number of Fourier modes
es = 1.46^2; % substrate permittivity
eg = 2.08^2; % grating permittivity
ec = 1.46^2; % cover permittivity
 DE = calc_emn_lam(no,0.45,eg,ec);
Rsum=zeros(Nl,Nt);
        Tsum=zeros(Nl,Nt);
for i=1:Nl
    for j=1:Nt
        wl = wll(i);
        kx0 = kxx0(j)*(ec^0.5);
        wv = 2*pi/wl; % wavevector
        % dimensionless variables
        kg = wl/gp;
        kh = wv*gh;
        
        
        
        %%
        
        
        % calculate in-plane wavevector projections and plane wave propagation
        % constants in the substrate (kz1) and in the cover (kz2)
        [kx,kz] = calc_kxz(no,kx0,kg,1.);
        [kz1,kz2] = calc_kzz(no,kx,es,ec);
        % calculate Fourier image of the dielectric permittivity function
       
        % scattering matrix of the grating
        SM = fmm(no,kh,kx,kz1,kz2,es,ec,DE,'TM');
        
        %%
        % diffraction of a plane wave
        VI = zeros(2,no); % first index indicate whether it is the substrate of the cover
        VI(2,floor(no/2)+1) = 1; % plane wave (0-th harmonic) coming from the cover
        VD = zeros(2,no);
        VD(1,:) = VI(1,:)*SM(:,:,1,1) + VI(2,:)*SM(:,:,2,1); % diffraction to the substrate
        VD(2,:) = VI(1,:)*SM(:,:,1,2) + VI(2,:)*SM(:,:,2,2); % diffraction to the cover
        % check energy conservation
        b = calc_balance(no,VI,VD,kz1,kz2,es,ec,'TM');
        disp(b);
        VE = calc_eff(no,VI,VD,kz1,kz2,es,ec,'TM');
        Tsum(i,j)=VE(1,floor(no/2)+1); % zero order power transmission
        Rsum(i,j)=VE(2,floor(no/2)+1); % zero order power reflection
    end
end
XI = wll;
YI = theta*180/pi;
ZI = transpose(Rsum);
figure;
pcolor(XI,YI,ZI)
shading flat
colormap('jet');
colorbar;
%caxis([0 0.4])
hold off
%%
