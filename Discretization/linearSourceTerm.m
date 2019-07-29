function M = linearSourceTerm(k)
% Matrix of coefficients for a linear source term in the form of k \phi
%
% k is a cell variable
%
% SYNOPSIS:
%   M = linearSourceTerm(k)
%
% PARAMETERS:
%
%
% RETURNS:
%
%
% EXAMPLE:
%
% SEE ALSO:
%

% Copyright (c) 2012-2016 Ali Akbar Eftekhari
% See the license file

d = k.domain.dimension;
if (d ==1) || (d==1.5) || (d==1.5)
	M = linearSourceTerm1D(k);
elseif (d == 2) || (d == 2.5) || (d==2.8)
	M = linearSourceTerm2D(k);
elseif (d == 3) || (d==3.2)
    M = linearSourceTerm3D(k);
end
