%% Propgagaion constant and mode field of GI-MM-Fibre
% This is a usage example of the mode solver. A graded-index multimmode
% fiber is considered in this example, which has been manufactured with a
% two-step thermal ion-exchange process using a borosilicate glass
% substrate and a silver nitrade melt solution. The fibure is highly
% multimodal. Any number of propagation constants can be obtained depending
% on the used computing hardware.

clear all
close all
clc

% Loading GI-Fiber data
load('diffused_waveguide.mat', 'nc');
n = nc;
n1 = max(max(n));  % Index of core 
n2 = min(min(n)); % Index of cladding
lambda = 1330e-9;       % Wavelength
beta_0 = 2*pi/lambda;   % Wave number
NA = sqrt(n1^2-n2^2);   % Maximum local numerical aperture

% Grid
x = -40:.1:40;
y = -10:.1:10;
[xg,yg] = meshgrid(x,y);
xg = xg*1e-6;
yg = yg*1e-6;

% Additional Parameters
POLARIZATION = 'TE';
FIELDCOMPONENT = 'Ex';
BC = 'ABC';
nbEigenwerte = 10;
nbInterpolations = 0;

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
eigenvalues = eigenvalues(1:ind(end));
n_eff = n_eff(1:ind(end));
modeFields = modeFields(:,1:ind(end));

out = ['The following effective indices have been calculated:'];
disp(out)
disp(num2str(n_eff,'%1.8f'));

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
surf(xgi,ygi,abs(modeFieldFundamental)/max(max(abs(modeFieldFundamental))))
shading interp
xlabel('x [um]')
ylabel('y [um]')
title('Normalized mode field of fundamental mode [a.u.]')