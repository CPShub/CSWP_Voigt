function [yi] = crossectional_stretch_sensitivity_entry(geo, mesh, mat, eps0, k0, u, K, iq)
% Compute the sensitivity of cross-sectional deformation u in relation to the
% strain prescriptor index (iq)
%
% Equations reference citation (1) as well as (3): "Numerische Methoden zur 
% Modellierung elastoplastischer Balken und ihre Anwendung auf 
% periodische Gitterstrukturen", PhD Thesis by L. Herrnböck
%
% Input:
	% geo   - Employed IGA Geometry 
	% mesh  - Employed mesh 
	% mat   - (Struct) containing material parameters
    % eps0  - (3,1) vector containing the strain prescriptors
    % k0    - (3,1) vector containing the twist prescriptors
    % u     - Displacement solution vector
    % K     - Stiffness matrix
    % iq    - Strain prescriptor index
% Output:
	% yi    - (n, 1) cross-sectional stretch sensitivity vector
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


% Implementation following Herrnboeck Equation 4.46

dof = 3;
ndofs = dof * mesh.nCpts;      % total dofs

gp_x = mesh.p+1;        % number of integration points in x-direction
gp_y = mesh.q+1;        % number of integration points in y-direction
[gp, wgt] = gauss_quadrature(gp_x, gp_y);   % calculate integration points and their weights

% Preallocate space
R_y = zeros(ndofs, 1);

% Compute derivatives of eps0and k0 -> selection vectors
[deps0_dq, dk0_dq] = derivatives_vk(eps0, k0, iq);

% Prenotate the unit tensor
e = eye(3);


% Assemlby of the R_y vector
count = 0;
for el = 1:mesh.nElems                % loop over elements
    sctr = mesh.elNodeCnt(el,:);       % element control points index
    elDoma = mesh.elDoma(el,:);        % element parametric domain
    elCpts0 = mesh.initcoords(sctr,:); % initial coordinates of el cont points
    nn = length(sctr);                % number of control points for each element
    nnElem = nn*dof;                  % dof for each element
    sctrB = zeros(1, nnElem);
    for i = 1:dof
        sctrB(i:dof:nnElem) = dof*(sctr-1) + i;  % displacement in i-th direction
    end
    
    elDisp = u(sctrB);
    elDisp = reshape(elDisp, dof, nn);
    elCpts(:,1:3)=elCpts0(:,1:3)+elDisp'; %here we actualize x+du
    
    for ipt = 1:size(gp,1)            % loop over integration points
        count = count + 1;
        pt = gp(ipt,:);      % reference parametric coordinates for each integration point
        wt = wgt(ipt);       % weigths for each integration point
        gauPts = parameter_gauss_mapping( elDoma, pt );   % gauss integration mapping
        j1 = jacobian_gauss_mapping( elDoma );     % jacobian value for gauss mapping   
        [N,ders] = nurbs_derivatives( gauPts,geo, mesh );
        jmatrix = ders*elCpts0(:,1:dof-1); %Because the mapping is in 2D
        j2 = det(jmatrix);
        ders =  jmatrix \ ders;              
        fac = j1 *j2 * wt;     
        ders3D = zeros(3,size(elCpts,1));
        ders3D(1:2,:) = ders;
        x = N.*elCpts(:,1:dof)'; % Position of the current Gauss Point in the mesh
        x = sum(x,2);
        x0 = N.*elCpts0(:,1:dof)';
        x0 = sum(x0,2);

        dx_alpha = elDisp * ders3D';
        F = def_gradient(eps0, k0, x, dx_alpha);

        if (det(F))<0
            warning('det(F) is negative, error')
        end

        if (mat.index >= 10 && mat.index < 20)
            % PK1-Based formulation
            [ pk1, A ] = material_CSWP_hyperelasticity(dof, mat, F); 
            
            % Equation 4.46 in (3)
            term1 = tensorprod(A, (deps0_dq + cross(dk0_dq, x)) * e(:,3)', [3 4], [1 2]);
            
            R_y_partial = zeros(nnElem, 1);
    
            countI = 1;
            for I=1:nn %loop over all NI --> Ry_I
                IDN = ders3D(:,I)';
                IN = N(I);
                for q=1:3
                    NI = IN*e(:,q);
                    DNI = zeros(3);
                    DNI(q,:)=IDN;
                    
                    term2 = cross(dk0_dq, pk1 * e(:,3))' * NI;
    
                    % Workaround because cross with term1 (3 3) and k0 (3 1)
                    % does not inherently work
                    term3_mat = zeros(3,3);
                    for j = 1:3
                        term3_mat(:,j) = cross(k0, term1(:,j));
                    end
                    term3 = (term3_mat * e(:, 3))' * NI;
    
                    R_y_partial(countI) = tensorprod(term1, DNI, "all") - term2 - term3;
                    countI = countI + 1;
                end
            end
        elseif (mat.index >= 110 && mat.index < 120)
            % PK2-Formulation, as presented in (1)
            [ pk2, dtan ] = material_CSWP_PK2_hyperelasticity(dof, mat, F);

            % B_u^{I} for PK2
            BN = zeros(6, nn*3);

            % B_uq^{I} for PK2
            BNq = zeros(nn*3, 6);

            term_a = (deps0_dq + cross(dk0_dq, x));
            dF_dq = term_a * e(3, :);
            EQ = to_voigt(0.5 * (F' * dF_dq + dF_dq' * F), "strain");
 
            for i = 1:nn  %loop over all NI --> Ry_I
                % See Eq. 85 in (1)
                col_3 = -N(i) * (cross(dk0_dq, F(:, 3)) + cross(k0, term_a));
                col_5 = ders(2, i) * term_a - N(i) * cross(dk0_dq, F(:, 2));
                col_6 = ders(1, i) * term_a - N(i) * cross(dk0_dq, F(:, 1));
                BNq(i*3-2:i*3, :) = [zeros(3, 1),zeros(3, 1),col_3,zeros(3, 1),col_5,col_6];

                % Corresponds to Eq. 52 in (1)
                BN(:,i*3-2:i*3) = [ F(1,1)*ders(1,i)     F(2,1)*ders(1,i)      F(3,1)*ders(1,i);
                    F(1,2)*ders(2,i)     F(2,2)*ders(2,i)      F(3,2)*ders(2,i);
                    N(i)*(k0(3)*F(2,3)-k0(2)*F(3,3))    N(i)*(k0(1)*F(3,3)-k0(3)*F(1,3))       N(i)*(k0(2)*F(1,3)-k0(1)*F(2,3)) ;
                    F(1,1)*ders(2,i)+ F(1,2)*ders(1,i)  F(2,1)*ders(2,i)+F(2,2)*ders(1,i)   F(3,1)*ders(2,i)+F(3,2)*ders(1,i);
                    (F(1,3)*ders(2,i) + N(i)*(k0(3)*F(2,2)-k0(2)*F(3,2)))  (F(2,3)*ders(2,i) + N(i)*(k0(1)*F(3,2)-k0(3)*F(1,2)))   (F(3,3)*ders(2,i) + N(i)*(k0(2)*F(1,2)-k0(1)*F(2,2)));
                    (F(1,3)*ders(1,i) + N(i)*(k0(3)*F(2,1)-k0(2)*F(3,1)))  (F(2,3)*ders(1,i) + N(i)*(k0(1)*F(3,1)-k0(3)*F(1,1)))   (F(3,3)*ders(1,i) + N(i)*(k0(2)*F(1,1)-k0(1)*F(2,1))) ];

            end
            % See Eq. 77 in (1)
            term1 = dtan * EQ;
            term2 = BNq * pk2;

            R_y_partial = (BN' * term1 + term2);
        end
        % Assembly into global R_y
        R_y(sctrB) = R_y(sctrB) + fac * R_y_partial;
    end
end

% Add 6 rows and columns of zeros as padding to R_y (see Eq 4.45 in (3))
R_y_padded = zeros(ndofs + 6, 1);
R_y_padded(1:ndofs, 1) = R_y;

% Solve for u_y
u_y_padded = -inv(K) * R_y_padded;

% Extract yi
yi = u_y_padded(1:ndofs);

end