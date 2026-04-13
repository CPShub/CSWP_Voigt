function [ pk2, dtan ] = constitutive_relation( dim, mat, F )
% Calculate constitutive relations with PK2
% Input:
    % dim   - (2 or 3) Dimension of the def. gradient
    % mat   - (Struct) containing the material parameters
    % F     - (2,2) or (3,3) def. gradient
% Output:
    % pk2   - (3,1) or (6,1) second piola-kirchhoff stress tensor in voigt
    % dtan  - (3,3) or (6,6) Material Elasticity Tensor
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

if (mat.index >= 10 && mat.index < 20)
    % belongs to hyperelastic materials with PK1 
    % -> Retrieve Model Response for index+100, as PK2 and dtan are
    % supposed to be provided
    mat2 = mat;
    mat2.index = mat.index + 100;
    
    if dim == 2
        mat2.compression_modulus = 0;
    end
    [ pk2, dtan ] = material_CWP_PK2_hyperelasticity( dim, mat2, F );
elseif (mat.index >= 110 && mat.index < 120)
    % belongs to hyperelastic materials with PK2
    mat2 = mat;
    if dim == 2
        %2d case is limited to incompressibility
        mat2.compression_modulus = 0;
    end
    [ pk2, dtan ] = material_CWP_PK2_hyperelasticity( dim, mat2, F );
elseif ( mat.index >= 20 && mat.index < 40 )        
    % belongs to plastic materials
elseif ( mat.index >= 40 && mat.index < 60 )        
    % belongs to viscous materials
end

end

