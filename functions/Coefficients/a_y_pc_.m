function [ a_y ] = a_y_pc_(n,yg,dim_xl,dim_yl,lok_adr_slgs,glob_adr_slgs,POLARIZATION,FIELDCOMPONENT)
    
% Calculate step size to 'norhtern' and 'southern' neighbor
n = yg(glob_adr_slgs - 1) - yg(glob_adr_slgs);
s = yg(glob_adr_slgs) - yg(glob_adr_slgs + 1);

% Round-off error
n = 1e-12*round(n*1e12);
s = 1e-12*round(s*1e12);

% Common coefficients
a_n = a_n_pc_(n,yg,dim_xl,dim_yl,lok_adr_slgs,glob_adr_slgs,POLARIZATION,FIELDCOMPONENT);
a_s = a_s_pc_(n,yg,dim_xl,dim_yl,lok_adr_slgs,glob_adr_slgs,POLARIZATION,FIELDCOMPONENT);

% Calculate coefficients
if (strcmp('TE',POLARIZATION)) == 1

    a_y = - a_n - a_s; 

elseif (strcmp('TM',POLARIZATION) && (strcmp('Hx',FIELDCOMPONENT) || strcmp('Hy',FIELDCOMPONENT))) == 1

    a_y = - a_n - a_s; 

elseif (strcmp('TM',POLARIZATION) && (strcmp('Ex',FIELDCOMPONENT) || strcmp('Ey',FIELDCOMPONENT))) == 1

    a_y = -4./(n.*s) + a_n + a_s;

else

    out = ['WARNING: Unknown field component or polarization. Possible choices are: ''Ex'', ''Ey'', ''Hx'', ''Hy'' or ''TE'' and ''TM'' respectively.'];
    disp(out);

end
   
end

