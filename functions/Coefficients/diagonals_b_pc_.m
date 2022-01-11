function [ Cb,Nb,Sb,Eb,Wb ] = diagonals_b_pc_( alpha,beta_0,beta_z,n_eff,ni,xg,yg,dz,dim_y,dim_xl,dim_yl,E_b,dGg,dGl,FIELDCOMPONENT,BPMTYPE,step)

% Calculating dielectric permittivity
epsilon = ni.*ni;

%% Defining address space
    
% Global element adresses in 'sparse SLE'
glob_adr_slgs   = dGg(2:end-1,2:end-1);
glob_adr_slgs   = reshape(glob_adr_slgs,size(glob_adr_slgs,1)*size(glob_adr_slgs,2),1);

% Global addresses of 'northern' neighbor
glob_adr_N = dGg(3:end-1,2:end-1); 
glob_adr_N = reshape(glob_adr_N,size(glob_adr_N,1)*size(glob_adr_N,2),1);

% Global addresses of 'southern' neighbor
glob_adr_S = dGg(2:end-2,2:end-1); 
glob_adr_S = reshape(glob_adr_S,size(glob_adr_S,1)*size(glob_adr_S,2),1);

% Global addresses of 'eastern' neighbor
glob_adr_E = dGg(2:end-1,2:end-2);
glob_adr_E = reshape(glob_adr_E,size(glob_adr_E,1)*size(glob_adr_E,2),1);

% Global addresses of 'western' neighbor
glob_adr_W = dGg(2:end-1,3:end-1);
glob_adr_W = reshape(glob_adr_W,size(glob_adr_W,1)*size(glob_adr_W,2),1);

% Local element adresses in 'sparse SLE'
lok_adr_slgs    = dGl(2:end-1,2:end-1);
lok_adr_slgs    = reshape(lok_adr_slgs,size(lok_adr_slgs,1)*size(lok_adr_slgs,2),1);

lok_adr_N  = dGl(3:end-1,2:end-1); 
lok_adr_N  = reshape(lok_adr_N,size(lok_adr_N,1)*size(lok_adr_N,2),1);

lok_adr_S  = dGl(2:end-2,2:end-1); 
lok_adr_S  = reshape(lok_adr_S,size(lok_adr_S,1)*size(lok_adr_S,2),1);

lok_adr_E = dGl(2:end-1,2:end-2);
lok_adr_E = reshape(lok_adr_E,size(lok_adr_E,1)*size(lok_adr_E,2),1);

lok_adr_W = dGl(2:end-1,3:end-1);
lok_adr_W = reshape(lok_adr_W,size(lok_adr_W,1)*size(lok_adr_W,2),1);

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

Cb =  E_b(glob_adr_slgs)          .* (a_x_(ni,xg,dim_xl,dim_yl,dim_y,lok_adr_slgs,glob_adr_slgs,POLARIZATION,FIELDCOMPONENT) + a_y_(ni,yg,dim_xl,dim_yl,lok_adr_slgs,glob_adr_slgs,POLARIZATION,FIELDCOMPONENT) + 2*1j*beta_z/(alpha*dz) + beta_0^2 .* (epsilon(glob_adr_slgs) - n_eff^2));
Nb =  E_b(glob_adr_slgs - 1)      .*  a_n_(ni,yg,dim_xl,dim_yl,lok_adr_N,glob_adr_N,POLARIZATION,FIELDCOMPONENT); 
Sb =  E_b(glob_adr_slgs + 1)      .*  a_s_(ni,yg,dim_xl,dim_yl,lok_adr_S,glob_adr_S,POLARIZATION,FIELDCOMPONENT);
Eb =  E_b(glob_adr_slgs + dim_y)  .*  a_e_(ni,xg,dim_y,dim_xl,dim_yl,lok_adr_E,glob_adr_E,POLARIZATION,FIELDCOMPONENT);
Wb =  E_b(glob_adr_slgs - dim_y)  .*  a_w_(ni,xg,dim_y,dim_xl,dim_yl,lok_adr_W,glob_adr_W,POLARIZATION,FIELDCOMPONENT);

end
