function [ a_x ] = a_x_pc_(n,xg,dim_xl,dim_yl,dim_y,lok_adr_slgs,glob_adr_slgs,POLARIZATION,FIELDCOMPONENT)

% Calculate step size to 'eastern' and 'western' neighbor
e = xg(glob_adr_slgs + dim_y) - xg(glob_adr_slgs);
w = xg(glob_adr_slgs) - xg(glob_adr_slgs - dim_y);

% Round-off error
e = 1e-12*round(e*1e12);
w = 1e-12*round(w*1e12);

% Common coefficients
a_w = a_w_pc_(n,xg,dim_y,dim_xl,dim_yl,lok_adr_slgs,glob_adr_slgs,POLARIZATION,FIELDCOMPONENT);
a_e = a_e_pc_(n,xg,dim_y,dim_xl,dim_yl,lok_adr_slgs,glob_adr_slgs,POLARIZATION,FIELDCOMPONENT);

% Calculate coefficients
if (strcmp('TM',POLARIZATION)) == 1

    a_x = - a_e - a_w; 

elseif (strcmp('TE',POLARIZATION) && (strcmp('Hx',FIELDCOMPONENT) || strcmp('Hy',FIELDCOMPONENT))) == 1

    a_x = - a_e - a_w; 

elseif (strcmp('TE',POLARIZATION) && (strcmp('Ex',FIELDCOMPONENT) || strcmp('Ey',FIELDCOMPONENT))) == 1

    a_x = -4./(e.*w) + a_e + a_w;

else

    out = ['WARNING: Unknown field component or polarization. Possible choices are: ''Ex'', ''Ey'', ''Hx'', ''Hy'' or ''TE'' and ''TM'' respectively.'];
    disp(out);

end
   
end

