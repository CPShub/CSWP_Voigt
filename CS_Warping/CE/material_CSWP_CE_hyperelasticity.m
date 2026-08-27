function [ pk2, dtan, dtan_c, dmu_dc, dmu_dC] = material_CSWP_CE_hyperelasticity(dim, mat, F, c, J_c, J_e, I1)
% General Framework for different chemo-elastic material models using PK2 Formulation
% Returns are in Voigt-Notation and follow engineering shear convention
% Input:
    % dim   - (2 or 3) Dimensionality of the mesh
    % mat   - (Struct) with index and relevant material properties
    % F     - (3,3) deformation gradient
    % c     - (Float) Local concentration
% Output:
    % For dim == 2
        % pk2   - (3,1) Second Piola–Kirchhoff stress (Voigt form)
        % dtan  - (3,3) Material Elasticity Tensor (Voigt form, engineering shear)
    % For dim == 3
        % pk2   - (6,1) Second Piola–Kirchhoff stress (Voigt form)
        % dtan  - (6,6) Material Elasticity Tensor (Voigt form, engineering shear)
% Info:
% Index | Name          | Energy-Formulations for the compressible 3d models
%   110 | Neo-Hook      | W(J1,J3) = A10*(J1-3)+K/2*(J3-1)^2
%   111 | Mooney-Rivlin | W(J1,J2,J3) = A10*(J1-3)+A01*(J2-3)+K/2*(J3-1)^2
%   112 | Yeoh          | W(J1,J3) = A10*(J1-3)+A20*(J1-3)^2+A30*(J1-3)^3+K/2*(J3-1)^2
%   114 | SVK (Pk2)     | - 
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

    e = eye(3,3);
    C = F'*F;
    iC = inv(C);
    J = det(F);
    J23 = J^(-2/3);
    ei1 = (e - 1/3 * I1 * iC);
    if dim == 2
        error("2D Chemoelastic material modeling not yet implemented");
    elseif dim == 3
        if mat.index == 310 % Chemo-hyperelastic Neo-Hookean
            K = mat.elastic_modulus1;
            G = mat.elastic_modulus2;
            
            pk2 = J_c * (K * J_e * (J_e - 1) * iC + G * J23 * ei1);

            % 2 * Partial derivative of Pk2 with respect to r. Cauchy-Green C
            % TODO: Double check that this is correct
            dtan = 2*G*J_c*J23 * ( ...
                (1/3)*I1*tensorProduct(iC,iC) ...
                - (1/3)*dyad(iC,e) ...
                - (1/3)*dyad(e,iC) ...
                + (1/9)*I1*dyad(iC,iC) ) ...
                + K*J*(2*J_e - 1)*dyad(iC,iC) ...
                - 2*K*J*(J_e - 1)*tensorProduct(iC,iC);
            
            % Partial derivative of Pk2 with respect to the concentration c
            dtan_c = G * mat.Omega * J23 * ei1 -K * mat.Omega * J_e^2 * iC;

            % Partial derivative of mu with respect to concentration c
            dmu_dc = (K * mat.Omega^2 * J^2) / (J_c^3);

            % Partial derivative of mu with respect to r. Cauchy-Green C
            dmu_dC = G * mat.Omega * J23 * ei1 - K * mat.Omega * J_e^2*iC;
        else
            error("No other chemo-hyperelastic material models currently available. Use Neo-Hookean (314)!")
        end
    end
end


