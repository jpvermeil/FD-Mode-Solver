function [ a_s_r ] = a_s_pc_(n,yg,dim_xl,dim_yl,lok_adr_S,glob_adr_S,POLARIZATION,FIELDCOMPONENT )

% Calculate step size to 'northern' and 'southern' neighbor
n = yg(glob_adr_S - 1) - yg(glob_adr_S);
s = yg(glob_adr_S) - yg(glob_adr_S + 1);

% Round-off error
n = 1e-12*round(n*1e12);
s = 1e-12*round(s*1e12);

% Calculating dielectric permittivity
epsilon = n.*n; 

% Calculate coefficients
if (strcmp('TE',POLARIZATION)) == 1

    a_s = (2./(s.*(s+n))) .* ones(length(glob_adr_S),1);

elseif (strcmp('TM',POLARIZATION) && (strcmp('Hx',FIELDCOMPONENT) || strcmp('Hy',FIELDCOMPONENT))) == 1

    a_s = (2./(s.*(s+n))) .* (2*epsilon(glob_adr_S) ./ (epsilon(glob_adr_S) + epsilon(glob_adr_S + 1))); 

elseif (strcmp('TM',POLARIZATION) && (strcmp('Ex',FIELDCOMPONENT) || strcmp('Ey',FIELDCOMPONENT))) == 1

    a_s = (2./(s.*(s+n))) .* (2*epsilon(glob_adr_S + 1) ./ (epsilon(glob_adr_S) + epsilon(glob_adr_S + 1)));

else

    out = ['WARNING: Unknown field component or polarization. Possible choices are: ''Ex'', ''Ey'', ''Hx'', ''Hy'' or ''TE'' and ''TM'' respectively.'];
    disp(out);

end

% Insert coefficients according to their correct address
a_s_r = zeros(dim_xl*dim_yl,1);
a_s_r(lok_adr_S) = a_s;
    
end