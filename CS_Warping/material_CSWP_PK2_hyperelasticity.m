function [pk2, dtan] = material_CSWP_PK2_hyperelasticity(dim, mat, F)
% General Framework for different material Models using PK2 Formulation
% Returns are in Voigt-Notation and follow engineering shear convention
% Input:
    % dim   - (2 or 3) Dimensionality of the mesh
    % mat   - (Struct) with index and relevant material properties
    % F     - (3,3) deformation gradient
% Output:
    % For dim == 2
        % pk1   - (3,1) Second Piola–Kirchhoff stress (Voigt form)
        % A     - (3,3) Material Elasticity Tensor (Voigt form, engineering shear)
    % For dim == 3
        % pk1   - (6,1) Second Piola–Kirchhoff stress (Voigt form)
        % A     - (6,6) Material Elasticity Tensor (Voigt form, engineering shear)
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





C = F'*F;      % right Cauchy-Green deformation tensor
if dim == 2    % plane element
    C11 = C(1,1); C12 = C(1,2); C21 = C(2,1); C22 = C(2,2);
    C33 = 1 / ( C11*C22 - C12*C21 );   % satisfy det(C) = 1
    I1 = C11+C22+C33;
%     I2 = C11*C22 - C12*C21 + (C11+C22)*C33;
    I1C = [1-C33^2*C22; 1-C33^2*C11; C33^2*C12];
    I2C = [C22+C33-C33^2*(C11+C22)*C22;  C11+C33-C33^2*(C11+C22)*C11; -C12+C33^2*(C11+C22)*C12];

    I1CC = C33^2*[     2*C33*C22^2,  2*C33*C11*C22-1,    -2*C33*C22*C12;
                   2*C33*C11*C22-1,      2*C33*C11^2,    -2*C33*C11*C12;
                    -2*C33*C22*C12,   -2*C33*C11*C12,   2*C33*C12^2+0.5;];

    I2CC = [           2*C33^3*C22^2*(C11+C22)-2*C33^2*C22,   1-2*C33^2*(C11+C22)+2*C33^3*(C11+C22)*C11*C22,            C33^2*C12-2*(C11+C22)*C33^3*C12*C22; ...
             1-2*C33^2*(C11+C22)+2*C33^3*(C11+C22)*C11*C22,             2*C33^3*C11^2*(C11+C22)-2*C33^2*C11,            C33^2*C12-2*(C11+C22)*C33^3*C12*C11; ...               
                       C33^2*C12-2*(C11+C22)*C33^3*C12*C22,             C33^2*C12-2*(C11+C22)*C33^3*C12*C11,  2*(C11+C22)*C33^3*C12^2+C33^2*(C11+C22)/2-0.5;];
    K = mat.compression_modulus;               
    if mat.index == 110 && K == 0  % incompressible Neo-Hookean
        A10 = mat.NH_coef1;
        pk2 = 2*A10*I1C;
        dtan = 4*A10*I1CC;
    elseif mat.index == 111 && K == 0 % incompressible Mooney-Rivlin
        A10 = mat.MR_coef1;
        A01 = mat.MR_coef2;
        pk2 = 2*A10*I1C + 2*A01*I2C;
        dtan = 4*A10*I1CC + 4*A01*I2CC;
    elseif mat.index == 112 && K == 0 % incompressible Yeoh
        A10 = mat.YEOH_coef1;
        A20 = mat.YEOH_coef2;
        A30 = mat.YEOH_coef3;
        pk2 = 2*A10*I1C + 4*A20*(I1-3)*I1C + 6*A30*(I1-3)^2*I1C;
        dtan = 4*A10*I1CC + 8*A20*(I1C*I1C') + 8*A20*(I1-3)*I1CC +...
               24*A30*(I1-3)*(I1C*I1C') + 12*A30*(I1-3)^2*I1CC;  
    elseif mat.index == 113 && K == 0 % incompressible Bidtanerman
        A10 = mat.BIDTANERMAN_coef1;
        A01 = mat.BIDTANERMAN_coef2;
        A20 = mat.BIDTANERMAN_coef3;
        A30 = mat.BIDTANERMAN_coef4;
        pk2 = 2*A10*I1C + 2*A01*I2C + 4*A20*(I1-3)*I1C + 6*A30*(I1-3)^2*I1C;
        dtan = 4*A10*I1CC + 4*A01*I2CC + 8*A20*(I1C*I1C') + 8*A20*(I1-3)*I1CC +...
            24*A30*(I1-3)*(I1C*I1C') + 12*A30*(I1-3)^2*I1CC; 
    elseif mat.index == 114 % Saint-Venant-Kirchhoff with PK2 (2D)
        mu = mat.mue_mat;
        lambda = mat.lambda_mat;
        
        E_GL = 0.5 * (C - eye(2));
        trE = trace(E_GL);
        
        % Vector output for pk2
        pk2 = [lambda*trE + 2*mu*E_GL(1,1);   % index 1
               lambda*trE + 2*mu*E_GL(2,2);   % index 2
               2*mu*E_GL(1,2)];                % index 3
        % Matrix output for dtan
        dtan = [lambda + 2*mu, lambda,         0;
                lambda,         lambda + 2*mu, 0;
                0,              0,             mu];
    end

elseif dim == 3   % solid element, consider the nearly-compressible material
    % Right cauchy-green deformation tensor
    C = F'*F;  
    C11 = C(1,1); C22 = C(2,2); C33 = C(3,3); C12 = C(1,2); C23 = C(2,3); C31 = C(3,1);

    % Three invariants of deformation tesnor
    I1 = trace(C);   
    I2 = 0.5*( trace(C)*trace(C)-trace(C*C) );
    I3 = det(C);

    % The third modified invariant
    J3 = sqrt(I3);

    % Voigt representation of first-order derivatives of the three invariants
    I1C = [1;1;1;0;0;0];
    I2C = [C22+C33; C11+C33; C11+C22; -C12; -C23; -C31];
    I3C = [C22*C33-C23*C23; C11*C33-C31*C31; C11*C22-C12*C12; ...
           C23*C31-C12*C33; C12*C31-C11*C23; C12*C23-C22*C31];

    % Matrix representation of second-order derivatives of the three invariants
    I1CC = zeros(6,6);
    I2CC = [0, 1, 1, 0, 0, 0;  1, 0, 1, 0, 0, 0;  1, 1, 0, 0, 0, 0;     ...
            0, 0, 0, -0.5, 0, 0;  0, 0, 0, 0, -0.5, 0;  0, 0, 0, 0, 0, -0.5];
    I3CC = [0, C33, C22, 0, -C23, 0;  C33, 0, C11, 0, 0, -C31;  C22, C11, 0, -C12, 0, 0;  ...
            0, 0, -C12, -C33/2, C31/2, C23/2;  -C23, 0, 0, C31/2, -C11/2, C12/2;  0, -C31, 0, C23/2, C12/2, -C22/2];   

    % First derivatives of the modified invariants   
    J1C = I3^(-1/3)*I1C - 1/3*I1*I3^(-4/3)*I3C;
    J2C = I3^(-2/3)*I2C - 2/3*I2*I3^(-5/3)*I3C;
    J3C = 0.5*I3^(-1/2)*I3C;  

    % Second-order derivatives of the modified invariants
    J1CC = I3^(-1/3)*I1CC + 4/9*I1*I3^(-7/3)*(I3C*I3C') -  1/3*I3^(-4/3)*( I1C*I3C'+I3C*I1C'+I1*I3CC );
    J2CC = I3^(-2/3)*I2CC + 10/9*I2*I3^(-8/3)*(I3C*I3C') - 2/3*I3^(-5/3)*( I2C*I3C'+I3C*I2C'+I2*I3CC );
    J3CC = 0.5*I3^(-0.5)*I3CC - 0.25*I3^(-3/2)*(I3C*I3C');

    K = mat.compression_modulus;
    if mat.index == 110  % Neo-Hookean
        % W(J1,J3) = A10*(J1-3)+K/2*(J3-1)^2
        A10 = mat.NH_coef1;
        pk2 = 2*( A10*J1C + K*(J3-1)*J3C );
        dtan = 4*(A10*J1CC + K*(J3C*J3C') + K*(J3-1)*J3CC);  
    elseif mat.index == 111 % Mooney-Rivlin
        %  W(J1,J2,J3) = A10*(J1-3)+A01*(J2-3)+K/2*(J3-1)^2
        A10 = mat.MR_coef1;
        A01 = mat.MR_coef2;
        pk2 = 2*( A10*J1C + A01*J2C + K*(J3-1)*J3C );
        dtan = 4*(A10*J1CC + A01*J2CC + K*(J3C*J3C') + K*(J3-1)*J3CC);        
    elseif mat.index == 112 % Yeoh
        % W(J1,J3) = A10*(J1-3)+A20*(J1-3)^2+A30*(J1-3)^3+K/2*(J3-1)^2
        A10 = mat.YEOH_coef1;
        A20 = mat.YEOH_coef2;
        A30 = mat.YEOH_coef3;
        pk2 = 2* ( A10*J1C +2*A20*(J1-3)*J1C+3*A30*(J1-3)^2*J1C+ K*(J3-1)*J3C );
        dtan = 2* ( A10*J1CC +2*A20*J1C*J1C+2*A20*(J1-3)*J1C+6*A30*(J1-3)*J1C*J1C + 3*A30(J1-3)^2*J1CC+ K*(J3C*J3C') + K*(J3-1)*J3CC );
    elseif mat.index == 113 % Bidtanerman
        % to be updated
    elseif mat.index == 114 % Saint-Venant-Kirchhoff with PK2 (3D)
        mu = mat.mue_mat;
        lambda = mat.lambda_mat;
        
        E_GL = 0.5 * (C - eye(3));
        trE = trace(E_GL);
        
        % Voigt Index Mapping: 11->1, 22->2, 33->3, 12->4, 23->5, 31->6
        pk2 = [lambda*trE + 2*mu*E_GL(1,1);   % index 1
               lambda*trE + 2*mu*E_GL(2,2);   % index 2
               lambda*trE + 2*mu*E_GL(3,3);   % index 3
               2*mu*E_GL(1,2);                % index 4
               2*mu*E_GL(2,3);                % index 5
               2*mu*E_GL(3,1)];                % index 6
               
        dtan = zeros(6,6);
        dtan(1:3, 1:3) = lambda + 2*mu*eye(3);
        dtan(1,2)=lambda; 
        dtan(1,3)=lambda; 
        dtan(2,1)=lambda; 
        dtan(2,3)=lambda; 
        dtan(3,1)=lambda; 
        dtan(3,2)=lambda;
        dtan(4,4)=mu; 
        dtan(5,5)=mu; 
        dtan(6,6)=mu;
    end
end

end


