function [n0, m0] = beam_forces(geo, mesh, mat, eps0, k0, u)
% Compute the three beam forces and moments acting on the cross-section
% Input:
	% geo   - Employed IGA Geometry 
	% mesh  - Employed mesh 
	% mat   - (Struct) containing material parameters
    % eps0  - (3,1) vector containing the strain prescriptors
    % k0    - (3,1) vector containing the twist prescriptors
    % u     - Displacement solution vector
% Output:
	% n0    - (3,1) forces acting on the cross-section (in kN)
	% m0    - (3,1) moments acting on the cross-section (in kNmm)
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
% (1) J.C. Alzate Cobo, T. Henkels and O. Weeger, "The cross-sectional 
% warping problem for hyperelastic beams: An efficient formulation in 
% Voigt notation", DOI: 10.48550/arXiv.2604.12886 
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


n0 = zeros(1, 3);
m0 = zeros(1, 3);

% Solve beam forces and moments for each entry in eps0(-> n0) and k0 (-> m0)
for i = 1:6
    result = beam_forces_entry(geo, mesh, mat, eps0, k0, u, i);
    if i < 4
        n0(i) = result;
    else
        m0(i-3) = result;
    end
end