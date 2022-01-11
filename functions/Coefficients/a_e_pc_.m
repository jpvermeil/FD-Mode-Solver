function [ a_e_r ] = a_e_pc_(n,xg,dim_y,dim_xl,dim_yl,lok_adr_E,glob_adr_E,POLARIZATION,FIELDCOMPONENT)

% Calculate step size to 'eastern' and 'western' neighbor
e = xg(glob_adr_E + dim_y) - xg(glob_adr_E);
w = xg(glob_adr_E) - xg(glob_adr_E - dim_y);

% Round-off error
e = 1e-12*round(e*1e12);
w = 1e-12*round(w*1e12);

% Calculating dielectric permittivity
epsilon = n.*n; 

% Calculate coefficients
if (strcmp('TM',POLARIZATION)) == 1

    a_e = 2./(e.*(e+w)); 

elseif (strcmp('TE',POLARIZATION) && (strcmp('Hx',FIELDCOMPONENT) || strcmp('Hy',FIELDCOMPONENT))) == 1

    a_e = (2./(e.*(e+w))) .* (2*epsilon(glob_adr_E) ./ (epsilon(glob_adr_E) + epsilon(glob_adr_E + dim_y))); 

elseif (strcmp('TE',POLARIZATION) && (strcmp('Ex',FIELDCOMPONENT) || strcmp('Ey',FIELDCOMPONENT))) == 1

    a_e = (2./(e.*(e+w))) .* (2*epsilon(glob_adr_E + dim_y) ./ (epsilon(glob_adr_E) + epsilon(glob_adr_E + dim_y)));

else

    out = ['WARNING: Unknown field component or polarization. Possible choices are: ''Ex'', ''Ey'', ''Hx'', ''Hy'' or ''TE'' and ''TM'' respectively.'];
    disp(out);

end

% Insert coefficients according to their correct address
a_e_r = zeros(dim_xl*dim_yl,1);
a_e_r(lok_adr_E) = a_e;

end

