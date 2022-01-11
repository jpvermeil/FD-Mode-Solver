function [ a_n_r ] = a_n_pc_(n,yg,dim_xl,dim_yl,lok_adr_N,glob_adr_N,POLARIZATION,FIELDCOMPONENT)

% Calculate step size to 'northern' and 'southern' neighbor
n = yg(glob_adr_N - 1) - yg(glob_adr_N);
s = yg(glob_adr_N) - yg(glob_adr_N + 1);

% Round-off error
n = 1e-12*round(n*1e12);
s = 1e-12*round(s*1e12);

% Calculating dielectric permittivity
epsilon = n.*n; 

% Calculate coefficients
if (strcmp('TE',POLARIZATION)) == 1

    a_n = (2./(n.*(n+s))) .* ones(length(glob_adr_N),1);

elseif (strcmp('TM',POLARIZATION) && (strcmp('Hx',FIELDCOMPONENT) || strcmp('Hy',FIELDCOMPONENT))) == 1

    a_n = (2./(n.*(n+s))) .* (2*epsilon(glob_adr_N) ./ (epsilon(glob_adr_N) + epsilon(glob_adr_N - 1))); 

elseif (strcmp('TM',POLARIZATION) && (strcmp('Ex',FIELDCOMPONENT) || strcmp('Ey',FIELDCOMPONENT))) == 1

    a_n = (2./(n.*(n+s))) .* (2*epsilon(glob_adr_N - 1) ./ (epsilon(glob_adr_N) + epsilon(glob_adr_N - 1)));

else

    out = ['WARNING: Unknown field component or polarization. Possible choices are: ''Ex'', ''Ey'', ''Hx'', ''Hy'' or ''TE'' and ''TM'' respectively.'];
    disp(out);

end

% Insert coefficients according to their correct address
a_n_r = zeros(dim_xl*dim_yl,1);
a_n_r(lok_adr_N) = a_n;

end

