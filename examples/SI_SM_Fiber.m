%% Propgagaion constant and mode field of SI-SM-Fibre
% This is a usage example of the mode solver. A step-index single-mode
% fiber is considered in this example. The propagation constant is compared
% to the analytical solution of the problem, which can be obtained using
% the normalized frequency for the given fiber specs.

clear all
close all
clc

% Definition of SI-Fiber
n1 = 1.45;  % Index of core 
n2 = 1.448; % Index of cladding
r = 5e-6;   % Radius
lambda = 1330e-9;       % Wavelength
beta_0 = 2*pi/lambda;   % Wave number
NA = sqrt(n1^2-n2^2);   % Numerical aperture
R = 2*pi*r*NA/lambda;   % Fiber parameter

% Grid
x = -15:.1:15;
y = -15:.1:15;
[xg,yg] = meshgrid(x,y);
xg = xg*1e-6;
yg = yg*1e-6;

% Refractive index profile
n = n2*ones(length(x),length(y));
n(sqrt(xg.^2+yg.^2) < r) = n1;

POLARIZATION = 'TE';
FIELDCOMPONENT = 'Ex';
BC = 'ABC';
nbEigenwerte = 3;
nbInterpolations = 0;

%% Analytical Solution

% Normalized Frequency (depending on Fiber parameter R!)
B = 0.3433;

beta_1 = beta_0 * n1;
beta_2 = beta_0 * n2;

beta_z = sqrt(B * (beta_1^2 - beta_2^2) + beta_2^2);
n_eff_analytical = beta_z/beta_0;

out = ['Analytical effective index of fundamental mode: ' num2str(n_eff_analytical,'%1.8f')];
disp(out)

%% Interpolation of index profile
% While the interpolation of the index profile can potentially reduce the
% error when numerically computing the propagation constant, the step-index
% profile is not unproblematic when interpolated. Use in step-index case
% with caution or be sure to use an appropriate interpolation method. Thus,
% this feature is rather used with continuous index profiles such as
% graded-index waveguides.

if nbInterpolations > 0 && mod(nbInterpolations,1) == 0
    
    ni = interpn(n,nbInterpolations);
    xgi = interpn(xg,nbInterpolations);
    ygi = interpn(yg,nbInterpolations);
    
elseif nbInterpolations == 0
    
    ni = n;
    xgi = xg;
    ygi = yg;
    
else
    
    out = ['Invalid specification of ''nbInterpolations''! Has to be 0 or a positive integer value.'];
    disp(out)
    
end

%% Numerical Mode Solver

% Dimensions and grid
dim_y   = size(ni,1);
dim_x   = size(ni,2);
dim_yl   = dim_y - 2;
dim_xl   = dim_x - 2;
dGl = zeros(size(ni,1),size(ni,2));
dGg = zeros(size(ni,1),size(ni,2));
dGl(2:end-1,2:end-1) = reshape(1:1:dim_xl*dim_yl',dim_yl,dim_xl);
dGg(1:end) = 1:1:length(dGg(1:end));

% Execute Mode Solver
[eigenvalues,n_eff,modeFields] = FDPropagationconstantsSemivec(ni,beta_0,xgi,ygi,dim_y,dim_xl,dim_yl,dGg,dGl,POLARIZATION,FIELDCOMPONENT,nbEigenwerte);

% Finding guided modes
ind = find(n_eff >= n2);

% Omit non-guided modes
eigenvalues = eigenvalues(1:ind);
n_eff = n_eff(1:ind);
modeFields = modeFields(:,1:ind);

out = ['Numerical effective index of fundamental mode: ' num2str(n_eff,'%1.8f')];
disp(out)
out = ['Resulting absolute error: ' num2str(n_eff-n_eff_analytical,'%1.8d') '.' newline 'Resulting relative error [%]: ' num2str(100*(n_eff-n_eff_analytical)/n_eff_analytical,'%1.8d') '.'];
disp(out)

% Mode field generation of fundamental mode
modeFieldFundamental = zeros(dim_y,dim_x);
modeFieldFundamental(2:end-1,2:end-1) = reshape(modeFields(:,1),dim_yl,dim_xl);

%% Visualization

subplot(1,2,1)
surf(xg,yg,n)
shading interp
xlabel('x [um]')
ylabel('y [um]')
title('Index Profile')

subplot(1,2,2)
surf(xgi,ygi,modeFieldFundamental/max(max(modeFieldFundamental)))
shading interp
xlabel('x [um]')
ylabel('y [um]')
title('Normalized mode field of fundamental mode [a.u.]')