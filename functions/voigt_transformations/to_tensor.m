function T = to_tensor(voigt, type)
% Transforms a voigt tensor of type "type" into its respective matrix form
% "type" == "strain" assumes engineering-strain and thus halves
% off-diagonal values
% Input:
	% voigt     - (m,1) or (m,m) Voigt Tensor 
	% type      - (String) being "strain" or other
% Output:
	% T         - (n,n) or (n,n,n,n) Tensor
% ----------------------------------------
% Copyright (C) 2026 Tobias Henkels and Juan C. Alzate Cobo. 
% 
% This code is an extension and modification of the NLIGA framework 
% originally developed by Du et al. (2020). 
% 
% ------------------------------------------------------------------------ 
% CITATION: 
% If you use this code for your research, please cite: 
% 
% (1) J.C. Alzate Cobo, T. Henkels and O. Weeger, "Efficient formulation of 
% the cross-sectional warping problem of hyperelastic 3D beams in Voigt 
% notation", DOI: 10.48550/arXiv.2604.12886 
% (2) X. Du, G. Zhao, W. Wang, M. Guo, R. Zhang, J. Yang, "NLIGA: A MATLAB 
% framework for nonlinear isogeometric analysis", Computer Aided 
% Geometric Design, 80, 101869, 2020. 
% https://doi.org/10.1016/j.cagd.2020.101869 
% ------------------------------------------------------------------------ 
% LICENSE: 
% This function is free software: you can redistribute it and/or modify it 
% under the terms of the GNU General Public License as published by the 
% Free Software Foundation, either version 3 of the License, or (at your 
% option) any later version. (GPL-3.0-or-later) 
% 
% This program is distributed in the hope that it will be useful, but 
% WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
% or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License 
% for more details. 
% ------------------------------------------------------------------------ 
% CONTACT: 
% - Tobias Henkels (tobias.henkels@stud.tu-darmstadt.de) 
% - Juan C. Alzate Cobo (alzate@cps.tu-darmstadt.de) 
% Technische Universität Darmstadt, Germany 
% ------------------------------------------------------------------------


if size(voigt, 1) == 3
    % (3,1) -> (2,2) OR
    % (3,3) -> (2,2,2,2)
    dim = 2;
elseif size(voigt, 1) == 6
    % (6, 1) -> (3,3) OR
    % (6,6) -> (3,3,3,3)
    dim = 3;
end

if dim == 2
    map = [1,1; 2,2; 1,2];
    T = zeros(2,2);
else
    map = [1,1; 2,2; 3,3; 1,2; 2,3; 3,1];
    T = zeros(3,3);
end

% Vector to 2nd order tensor
if isvector(voigt)
    for a = 1:length(voigt)
        i = map(a,1); j = map(a,2);
        val = voigt(a);
        if strcmp(type, 'strain') && i ~= j
            T(i,j) = 0.5 * val;
            T(j,i) = 0.5 * val;
        else
            T(i,j) = val;
            T(j,i) = val;
        end
    end

% Matrix to 4th order tensor
else
    num_comp = size(voigt, 1);
    T = zeros(dim, dim, dim, dim);
    for a = 1:num_comp
        for b = 1:num_comp
            i = map(a,1); j = map(a,2);
            k = map(b,1); l = map(b,2);
            val = voigt(a,b);
            % Fill all 4 symmetric positions
            T(i,j,k,l) = val;
            T(j,i,k,l) = val;
            T(i,j,l,k) = val;
            T(j,i,l,k) = val;
        end
    end
end
end