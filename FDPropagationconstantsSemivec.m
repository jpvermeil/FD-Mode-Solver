function [ eigenvalues,n_eff,modeFields ] = FDPropagationconstantsSemivec(n,beta_0,xg,yg,dim_y,dim_xl,dim_yl,dGg,dGl,POLARIZATION,FIELDCOMPONENT,nbEigenwerte)
% Finite difference Mode Solver for TE/TM E- and/or H-fields in optical
% waveguide structures with arbitrary index profile. Calculation scheme is
% based on a semi-vectorial finite difference approach utilizing an
% absorbing boundary condition.
% 
% SYNOPSIS
%
% FDPropagationconstantsSemivec(n,beta_0,xg,yg,dim_y,dim_xl,dim_yl,dGg,dGl,POLARIZATION,FIELDCOMPONENT,NbEigenwerte)
% 
% VARIABLES
%
%   n               Refractive index profile.
%   beta_0          Wave vector according to specified wavelength.
%   xg              Grid of x dimension. Has to match dimensions of n. 
%                   Has to be X output of Matlab meshgrid function ([X,Y,Z] = meshgrid(x,y,z))
%   yg              Grid of y dimension. Has to match dimensions of n.
%                   Has to be Y output of Matlab meshgrid function ([X,Y,Z] = meshgrid(x,y,z))
%   dim_y           Global matrix dimension in second direction (see examples)
%   dim_xl          Local matrix dimension (without boundary values) in first direction (see examples)
%   dim_yl          Local matrix dimension (without boundary values) in second direction (see examples)
%   dGg             Global address vector grid
%   dGl             Local address vector grid
%   POLARIZATION    String value that can either be 'TE' or 'TM' (string)
%   FIELDCOMPONENTS Can be 'Ex', 'Ey', 'Hx' or 'Hy' (string)
%   nbEigenvalues   Desired number of eigenvalues that are to be determined. Caution: Not all of the calculated values are necessarily guided modes. Manual check required. 

% Finite difference Mode Solver for TE/TM E- and/or H-fields in optical
% waveguide structures with arbitrary index profile. Copyright (C) 2021
% Jan-Philipp Roth (jan-philipp.roth@protonmail.com) This program is free
% software; you can redistribute it and/or modify it under the terms of the
% GNU General Public License as published by the Free Software Foundation;
% either version 3 of the License, or (at your option) any later version.
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
% Public License for more details. You should have received a copy of the
% GNU General Public License along with this program; if not, write to the
% Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
% MA 02110-1301  USA

%% Initial Definitions, Checks and Parameters

if ((strcmp('TE',POLARIZATION) || strcmp('TM',POLARIZATION)) && ((strcmp('Ex',FIELDCOMPONENT) || strcmp('Ey',FIELDCOMPONENT) || strcmp('Hx',FIELDCOMPONENT) || strcmp('Hy',FIELDCOMPONENT)))) == 0

    out = 'WARNING: Unknown field component or polarization. Possible choices are: ''Ex'', ''Ey'', ''Hx'', ''Hy'' or ''TE'' and ''TM'' respectively.';
    disp(out);
    return

end

%% Generation of Diagonals

% The functions generating the coefficients a_*_ always yield vector that
% match the provided addresses. Hence, the format of addresses and
% coefficients does match. The adress vector contains the information for
% which element of the global system matrix the according coefficient is
% needed.

% The generation of the system matrix is only defined by the diagonal
% vectors and not by the characteristic of the coefficient itself. The
% coefficients only depend on polarization and field component and can be
% chosen.

% Calculation of the diagonal vectors
[ diagC,diagN,diagS,diagE,diagW ] = diagonals_x_pc_( beta_0,n,xg,yg,dim_y,dim_xl,dim_yl,dGg,dGl,POLARIZATION,FIELDCOMPONENT);

% Assembling of the system matrix with diagonal vectors
x = sparse((size(n,1)-2)*(size(n,2)-2),(size(n,1)-2)*(size(n,2)-2));
x = spdiags(diagC,0,x);
x = spdiags([diagN(2:end); 0],-1,x);
x = spdiags([0; diagS(1:end-1)],1,x);
x = spdiags([zeros(dim_yl,1); diagE(1:end-dim_yl)],dim_yl,x);
x = spdiags([diagW(dim_yl+1:end); zeros(dim_yl,1)],-dim_yl,x);

%% Calculation of the 'right side'

b = sparse(dim_yl*dim_xl,dim_yl*dim_xl);
b = spdiags(ones(dim_yl*dim_xl,1),0,b);

%% Solution of the Eigenvalue problem (xx - beta^2*b = 0)

opts.maxit = 1000;
opts.tol   = 1e-9;

[modeFields,eigenvalues] = eigs(x,b,nbEigenwerte,(max(max(n))*beta_0)^2,opts);

% Calculation of the effective indices n_eff
n_eff = sqrt(eigenvalues./beta_0^2);
n_eff = diag(n_eff,0);
n_eff = n_eff(1:nbEigenwerte);

eigenvalues = diag(eigenvalues);

% Calculation of mode fields (which correspond to the eigenvectors of
% the problem)
modeFields = modeFields(:,1:nbEigenwerte);
    
end

