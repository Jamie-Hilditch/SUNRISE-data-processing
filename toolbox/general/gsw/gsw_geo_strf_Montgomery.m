function geo_strf_Montgomery = gsw_geo_strf_Montgomery(SA,CT,p,p_ref)

% gsw_geo_strf_Montgomery                            Montgomery geostrophic 
%                                         streamfunction (75-term equation)
%==========================================================================
%
% USAGE:  
% geo_strf_Montgomery = gsw_geo_strf_Montgomery(SA,CT,p,p_ref)
%
% DESCRIPTION:
%  Calculates the Montgomery geostrophic streamfunction (see Eqn. (3.28.1) 
%  of IOC et al. (2010)).  This is the geostrophic streamfunction for the 
%  difference between the horizontal velocity at the pressure concerned, p,
%  and the horizontal velocity on the pressure surface, p_ref.  The 
%  Montgomery geostrophic streamfunction is the geostrophic streamfunction
%  for flow in a specifc volume anomaly surface.  The reference values used
%  for the specific volume anomaly are SA = SSO = 35.16504 g/kg and  
%  CT = 0 deg C.  This function calculates specific volume anomaly using 
%  the computationally efficient 75-term expression for specific volume of 
%  Roquet et al. (2015).
%
%  Note that p_ref, is the reference pressure to which the streamfunction
%  is referenced.  When p_ref is zero, "gsw_geo_strf_Montgomery" returns 
%  the Montgomery geostrophic streamfunction with respect to the sea 
%  surface, otherwise, the function returns the geostrophic streamfunction
%  with respect to the (deep) reference pressure p_ref.
%
%  Note that the 75-term equation has been fitted in a restricted range of 
%  parameter space, and is most accurate inside the "oceanographic funnel" 
%  described in McDougall et al. (2003).  The GSW library function 
%  "gsw_infunnel(SA,CT,p)" is avaialble to be used if one wants to test if 
%  some of one's data lies outside this "funnel".  
%
% INPUT:
%  SA    =  Absolute Salinity                                      [ g/kg ]
%  CT    =  Conservative Temperature (ITS-90)                     [ deg C ]
%  p     =  sea pressure                                           [ dbar ]
%           ( i.e. absolute pressure - 10.1325 dbar )
%  p_ref =  reference pressure                                     [ dbar ]
%           ( i.e. reference absolute pressure - 10.1325 dbar )
%
%  SA & CT need to have the same dimensions.
%  p may have dimensions Mx1 or 1xN or MxN, where SA & CT are MxN.
%  p_ref needs to be a single value, it can have dimensions 1x1 or Mx1 or  
%  1xN or MxN.
%
% OUTPUT:
%  geo_strf_Montgomery  =  Montgomery geostrophic               [ m^2/s^2 ]
%                          streamfunction                      
%
% AUTHOR:  
%  Trevor McDougall and Paul Barker                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.05 (27th January 2015)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See section 3.28 of this TEOS-10 Manual. 
%
%  McDougall, T.J., D.R. Jackett, D.G. Wright and R. Feistel, 2003: 
%   Accurate and computationally efficient algorithms for potential 
%   temperature and density of seawater.  J. Atmosph. Ocean. Tech., 20,
%   pp. 730-741.
%
%  Montgomery, R. B., 1937: A suggested method for representing gradient 
%   flow in isentropic surfaces.  Bull. Amer. Meteor. Soc. 18, 210-212.  
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 4)
    error('gsw_geo_strf_Montgomery: Requires four inputs')
end 

unique_p_ref = unique(p_ref);
if ~isscalar(unique_p_ref)
    error('gsw_geo_strf_Montgomery: The reference pressure p_ref must be unique')
end
clear p_ref
p_ref = unique_p_ref;

[ms,ns] = size(SA);
[mt,nt] = size(CT);
[mp,np] = size(p);

if (ms~=mt) | (ns~=nt)
    error('gsw_geo_strf_Montgomery: SA & CT need to have the same dimensions')
end

if (mp == 1) & (np == 1)              % p is a scalar 
    error('gsw_geo_strf_Montgomery: need more than one pressure'); 
elseif (ns == np) & (mp == 1)         % p is row vector,
    p = p(ones(1,ms), :);              % copy down each column.
elseif (ms == mp) & (np == 1)         % p is column vector,
    p = p(:,ones(1,ns));               % copy across each row.
elseif (ns == mp) & (np == 1)          % p is a transposed row vector,
    p = p.';                              % transposed then
    p = p(ones(1,ms), :);                % copy down each column.
elseif (ms == mp) & (ns == np)
    % ok
else
    error('gsw_geo_strf_Montgomery: Inputs array dimensions arguments do not agree')
end 

transposed = 0;
if ms == 1  
   p  =  p(:);
   CT  =  CT(:);
   SA  =  SA(:);
   transposed = 1;
end 

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

db2Pa = 1e4;

dyn_height  = gsw_geo_strf_dyn_height(SA,CT,p,p_ref);

geo_strf_Montgomery = db2Pa*p.*gsw_specvol_anom_standard(SA,CT,p) + dyn_height;
                             
%--------------------------------------------------------------------------
% This function calculates the Montgomery streamfunction using the 
% computationally efficient 75-term expression for specific volume in terms
% of SA, CT and p.  If one wanted to compute this with the full TEOS-10 
% Gibbs function expression for specific volume, the following lines of 
% code will enable this.  Note that dynamic height will also need to be 
% evaluated using the full Gibbs function.
%
%    geo_strf_Montgomery = db2Pa*p.*gsw_specvol_anom_standard_CT_exact(SA,CT,p) ...
%                               + dyn_height;
%
%---------------This is the end of the alternative code--------------------

if transposed
   geo_strf_Montgomery = geo_strf_Montgomery.';
end

end
