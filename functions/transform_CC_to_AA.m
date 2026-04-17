function [AA] = transform_CC_to_AA(F, S, CC)
% This function transforms the material elasticity tensor CC (also called
% dtan) into the First Elasticity tensor (AA).
% Input:
    % F     - (3,3) deformation gradient
    % S     - (3,3) Second Piloa-Kirchhoff Stress tensor
    % CC    - (3,3,3,3) Material Elasticity Tensor
% Output:
    % AA    - (3,3,3,3) First Elasticity Tensor

% Transformation follows equation:
%     frac{\partial PK1_{i,J}}{\partial F_{k,L}} 
%   = \mathbb{A}_{iJkL} 
%   = \delta_{ik} S_{JL} + F_{iM} \mathbb{C}_{MJNL} F_{kN}

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
    

AA = zeros(3,3,3,3);
delta = eye(3);

for J = 1:3
    for L = 1:3
        CC_JL = reshape(CC(:,J,:,L), [3,3]); 
        % Compute F_iM * CC_MJNL * F_kN for all i, k
        term_mat_ik = F * CC_JL * F'; 
        
        for i = 1:3
            for k = 1:3
                AA(i,J,k,L) = delta(i,k) * S(J,L) + term_mat_ik(i,k);
            end
        end
    end
end
end