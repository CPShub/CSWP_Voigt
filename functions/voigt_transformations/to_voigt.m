function voigt = to_voigt(T, type)
% Transforms a tensor of type "type" into its respective voigt form
% "type" == "strain" assumes engineering-strain and thus doubles
% off-diagonal values
% Input:
	% T      - (n,n) or (n,n,n,n) Tensor 
	% type   - (String) being "strain" or other
% Output:
	% voigt  - (m,1) or (m,m) Voigt Tensor
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
% notation", [Journal Name], [Year]. DOI: [DOI] 
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

% type: 'stress' or 'strain'
dim = size(T, 1);
tol = 1e-12;

% Symmetry check for 2nd order
if ndims(T) == 2
    if max(max(abs(T - T'))) > tol
        error('Tensor is not symmetric within tolerance');
    end
end

if dim == 2
    map = [1,1; 2,2; 1,2];
    num_comp = 3;
else
    map = [1,1; 2,2; 3,3; 1,2; 2,3; 3,1];
    num_comp = 6;
end

% 2nd order tensor to vector
if ndims(T) <= 2 
    voigt = zeros(num_comp, 1);
    for a = 1:num_comp
        i = map(a,1); j = map(a,2);
        val = T(i,j);
        if strcmp(type, 'strain') && i ~= j % Engineering-Strain assumes doubling of off-diagonals in voigt-notation
            voigt(a) = 2 * val;
        else
            voigt(a) = val;
        end
    end

% 4th order tensor to matrix
else 
    voigt = zeros(num_comp, num_comp);
    for a = 1:num_comp
        for b = 1:num_comp
            i = map(a,1); j = map(a,2);
            k = map(b,1); l = map(b,2);
            % CC is copied 1:1 regardless of type 
            % if strains are scaled by 2 in the strain-vector
            voigt(a,b) = T(i,j,k,l);
        end
    end
end
end