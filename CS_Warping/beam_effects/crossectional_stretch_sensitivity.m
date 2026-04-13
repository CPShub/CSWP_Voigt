function [y] = crossectional_stretch_sensitivity(geo, mesh, mat, eps0, k0, u, K)
% Compute the sensitivity of cross-sectional deformation u in relation to the
% 6 strain prescriptors
%
% Equations reference "Numerische Methoden zur 
% Modellierung elastoplastischer Balken und ihre Anwendung auf 
% periodische Gitterstrukturen", PhD Thesis by L. Herrnböck
%
% Input:
	% geo   - Employed IGA Geometry 
	% mesh  - Employed mesh 
	% mat   - struct containing material parameters
    % eps0  - (3,1) vector containing the strain prescriptors
    % k0    - (3,1) vector containing the twist prescriptors
    % u     - displacement solution vector
    % K     - stiffness matrix
% Output:
	% y    - (n,6) cross-sectional stretch sensitivity matrix
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

% Solver for crossectional stretch senstivity following
% Herrnboeck Equations 4.43 to 4.46

dof = 3;
ndofs = dof * mesh.nCpts;      % total dofs

% Preallocate space;
y = zeros(ndofs, 6);

% Solve for crossectional stretch sensitivity for each entry q in p
for iq = 1:6 % parfor?
    [yi] = crossectional_stretch_sensitivity_entry(geo, mesh, mat, eps0, k0, u, K, iq);

    y(:, iq) = yi;
end

end
