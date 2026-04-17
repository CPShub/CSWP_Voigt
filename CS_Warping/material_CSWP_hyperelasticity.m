function [pk1, A] = material_CSWP_hyperelasticity (dim, mat, F)
% General Framework for different material models using PK1 Formulation
% Input:
    % dim   - (2 or 3) Dimensionality of the mesh
    % mat   - (Struct) with index and relevant material properties
    % F     - (3,3) deformation gradient
% Output:
    % For dim == 2
        % pk1   - (2,2) First Piola–Kirchhoff stress
        % A     - (2,2,2,2) First Elasticity Tensor
    % For dim == 3
        % pk1   - (3,3) First Piola–Kirchhoff stress
        % A     - (3,3,3,3) First Elasticity Tensor
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





% Index-Associations (mat.index)
%   10: Neo-Hooke
%   11: Mooney-Rivli
%   12: Yeoh
%   13: Bidtanerman
%   14: Saint-Venant-Kirchhoff (PK1)

% Change material index to line up with PK2 Formulations
mat2 = mat;
mat2.index = mat.index + 100;

if dim == 2
    % Start by Computing pk2 and Material Elasticity Tensor dtan
    [pk2, dtan] = material_CSWP_PK2_hyperelasticity(dim, mat2, F);
    
    pk2_full = to_tensor(pk2, "stress");
    dtan_full = to_tensor(dtan, "-");

    % Shift from PK2 and dtan=2*Derivative(PK2, C)=Derivative(PK2, E) 
    % to PK1 and A=Derivative(PK1, F)
    term_mat_full = tensorprod(F, dtan_full, 2, 1);     
    term_mat_full = tensorprod(term_mat_full, F, 3, 2);  % Shifts Index Order
    % Reorder to (i,j,k,l)
    term_mat_full = permute(term_mat_full, [1, 2, 4, 3]);
    A = term_mat_full;
    for i = 1:2
        for j = 1:2
            for k = 1:2
                for l = 1:2 
                    % term_geo = delta_{ik} * S_{jl}
                    if i == k
                        A(i,j,k,l) = A(i,j,k,l) + pk2_full(j,l);
                    end
                end
            end
        end
    end
elseif dim == 3
    % Start by Computing pk2 and dtan
    [pk2, dtan] = material_CSWP_PK2_hyperelasticity(dim, mat2, F);
    
    pk2_full = to_tensor(pk2, "stress");
    dtan_full = to_tensor(dtan, "-");

    pk1 = F * pk2_full;
    A = transform_CC_to_AA(F, pk2_full, dtan_full);
end
















