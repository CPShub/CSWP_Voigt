function [C0] = beam_stiffness(geo, mesh, mat, eps0, k0, u, K)
% Compute the (6,6) beam stiffness matrix relating changed in strain
% prescriptors to changes in beam forces / moments
% Input:
	% geo   - Employed IGA Geometry 
	% mesh  - Employed mesh 
	% mat   - struct containing material parameters
    % eps0  - (3,1) vector containing the strain prescriptors
    % k0    - (3,1) vector containing the twist prescriptors
    % u     - displacement solution vector
    % K     - stiffness matrix
% Output:
	% C0    - (6,6) beam stiffness matrix
% ------------------------------------------------------------------------
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
C0 = zeros(6);

% Solve crossectional stretch sensitivity for each strain prescriptor
[y] = crossectional_stretch_sensitivity(geo, mesh, mat, eps0, k0, u, K);

for i = 1:6 % Counter for p (selection pointer for forces / moments)
    for j = 1:6 % Counter for q (selection pointer for strain prescriptors)
        C0(i, j) = beam_stiffness_entry(geo, mesh, mat, eps0, k0, u, i, j, y);
    end
end
end