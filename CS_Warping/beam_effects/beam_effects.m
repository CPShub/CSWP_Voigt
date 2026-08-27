function [n0, m0, C0, uy_all] = beam_effects(geo, mesh, mat, eps0, k0, u, K)
    % This function computes the beam forces, moments and beam stiffness, 
    % as well as the deformation solution sensitivities with regard to the 
    % applied strain prescriptors from the deformation solution.
    % Input:
        % geo           - Employed IGA Geometry
        % mesh          - Employed mesh 
        % mat           - (Struct) containing material parameters
        % eps0          - Vector containing the strain prescriptors
        % k0            - Vector containing the twist prescriptors
        % u             - Displacement solution vector
        % K             - Tangent Stiffness matrix for displacement solution vector u  

    % Output:
        % n0            - (3,1) Beam forces
        % m0            - (3,1) Beam moments
        % C0            - (6,6) Beam stiffnesses
        % uy_all        - Deformation solution sensitivities u,q
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

    
    % Integration over the whole domain
    gp_x = mesh.p+1;        % number of integration points in x-direction
    gp_y = mesh.q+1;        % number of integration points in y-direction
    [gp, wgt] = gauss_quadrature(gp_x, gp_y);   % calculate integration points and their weights
    
    dof = 3;
    ndofs = dof * mesh.nCpts;      % total dofs
    
    % Preallocate space
    R_y_all = zeros(ndofs, 6);
    stress_resultant = zeros(6,1);
   
    dk0dq = [zeros(3,3), eye(3,3)];
    deps0dq = [eye(3,3), zeros(3,3)];
    
    % First Integration to Compute the Stress Resultants n0, m0 and the Displacement
    % Sensitivites u,y
    for el = 1:mesh.nElems                % loop over elements
        sctr = mesh.elNodeCnt(el,:);       % element control points index
        elDoma = mesh.elDoma(el,:);        % element parametric domain
        elCpts0 = mesh.initcoords(sctr,:); % initial coordinates of el cont points
        nn = length(sctr);                % number of control points for each element
        nnElem = nn*dof;                  % total dof for each element
        sctrB = zeros(1, nnElem);
    
        for i = 1:dof
            sctrB(i:dof:nnElem) = dof*(sctr-1) + i;  % displacement in i-th direction
        end
        
        elDisp = u(sctrB);
        elDisp = reshape(elDisp, dof, nn);
        elCpts(:,1:3)=elCpts0(:,1:3)+elDisp'; %here we actualize x+du
        
        for ipt = 1:size(gp,1)            % loop over integration points
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
    
            dx_alpha = elDisp * ders3D';
            F = def_gradient(eps0, k0, x, dx_alpha);
    
            if (det(F))<0 % Check determinate 
                warning('det(F) is negative, error')
            end
            
            % Compute the material response
            if (mat.index >= 10 && mat.index < 20) % Formulation with PK1
                error("This solution scheme for the beam effects and stiffness is only usable with a PK2 material formulation!")
            elseif (mat.index >= 110 && mat.index < 120) % Formulation with PK2
                [pk2, dtan] = material_CSWP_PK2_hyperelasticity(dof, mat, F);
            end
    
            % % Compute BN
            % BN = zeros(6, nn*3);
            % 
            % % Indicees for blocks
            % idx1 = 1:3:3*nn;
            % idx2 = 2:3:3*nn;
            % idx3 = 3:3:3*nn;
            % 
            % C3 = [ (k0(3)*F(2,3)-k0(2)*F(3,3)), (k0(1)*F(3,3)-k0(3)*F(1,3)), (k0(2)*F(1,3)-k0(1)*F(2,3)) ];
            % C5 = [ (k0(3)*F(2,2)-k0(2)*F(3,2)), (k0(1)*F(3,2)-k0(3)*F(1,2)), (k0(2)*F(1,2)-k0(1)*F(2,2)) ];
            % C6 = [ (k0(3)*F(2,1)-k0(2)*F(3,1)), (k0(1)*F(3,1)-k0(3)*F(1,1)), (k0(2)*F(1,1)-k0(1)*F(2,1)) ];
            % 
            % 
            % % Row 1
            % BN(1, idx1) = F(1,1)*ders(1,:); 
            % BN(1, idx2) = F(2,1)*ders(1,:); 
            % BN(1, idx3) = F(3,1)*ders(1,:);
            % 
            % % Row 2
            % BN(2, idx1) = F(1,2)*ders(2,:); 
            % BN(2, idx2) = F(2,2)*ders(2,:); 
            % BN(2, idx3) = F(3,2)*ders(2,:);
            % 
            % % Row 3
            % BN(3, idx1) = C3(1)*N;          
            % BN(3, idx2) = C3(2)*N;          
            % BN(3, idx3) = C3(3)*N;
            % 
            % % Row 4
            % BN(4, idx1) = F(1,1)*ders(2,:) + F(1,2)*ders(1,:);
            % BN(4, idx2) = F(2,1)*ders(2,:) + F(2,2)*ders(1,:);
            % BN(4, idx3) = F(3,1)*ders(2,:) + F(3,2)*ders(1,:);
            % 
            % % Row 5
            % BN(5, idx1) = C5(1)*N + F(1,3)*ders(2,:);
            % BN(5, idx2) = C5(2)*N + F(2,3)*ders(2,:);
            % BN(5, idx3) = C5(3)*N + F(3,3)*ders(2,:);
            % 
            % % Row 6
            % BN(6, idx1) = F(1,3)*ders(1,:) + C6(1)*N;
            % BN(6, idx2) = F(2,3)*ders(1,:) + C6(2)*N;
            % BN(6, idx3) = F(3,3)*ders(1,:) + C6(3)*N;
            % 
    
            % Compute Derivatives of BN wrt q
            R_y_partial_all = zeros(nnElem, 6);
            
            % Extract the relevant stress components for PK2_reduced
            S_red = [pk2(6); pk2(5); pk2(3)];
    
            % Extract the relevant columns of the material stiffness matrix
            C_red = dtan(:, [6, 5, 3]); %dtan(:, [3,5,6]);
    
            % Compute explicit strain derivatives
            dkdq_cross_x = [0,x(3),-x(2); -x(3), 0, x(1); x(2), -x(1), 0];
            deps0k0_dq = [eye(3,3), dkdq_cross_x];

            % compute the h,p operators for all variations of p in [eps0, k0]
            % Assemble them directly into the H operator
            H_eps0 = F';
            H_k0 = F' * dkdq_cross_x;
            H = [H_eps0, H_k0];

            % Compute Combinations of strain prescriptor derivatives
            dkdq_cross_F1 = [zeros(3,3), cross(eye(3,3), repmat(F(:,1), 1, 3))];
            dkdq_cross_F2 = [zeros(3,3), cross(eye(3,3), repmat(F(:,2), 1, 3))];
            dkdq_cross_F3 = [zeros(3,3), cross(eye(3,3), repmat(F(:,3), 1, 3))];
   
            % Compute the E,q^{0} operator following eq. 118 in (19
            E_v_0_q = zeros(6,6);

            % Only fill the 3rd, 5th and 6th layers
            E_v_0_q(3, :) = F(:, 3)' * deps0k0_dq;
            E_v_0_q(5, :) = F(:, 2)' * deps0k0_dq;
            E_v_0_q(6, :) = F(:, 1)' * deps0k0_dq;

            BN = zeros(6, nn*3);
            BN_6_q = zeros(nn*3, 6);
    
            for i = 1:nn
                % Compute BN for all N=1:nn
                BN(:,i*3-2:i*3) = [ 
                    F(1,1)*ders(1,i)     F(2,1)*ders(1,i)      F(3,1)*ders(1,i);
                    F(1,2)*ders(2,i)     F(2,2)*ders(2,i)      F(3,2)*ders(2,i);
                    N(i)*(k0(3)*F(2,3)-k0(2)*F(3,3))    N(i)*(k0(1)*F(3,3)-k0(3)*F(1,3))       N(i)*(k0(2)*F(1,3)-k0(1)*F(2,3)) ;
                    F(1,1)*ders(2,i)+ F(1,2)*ders(1,i)  F(2,1)*ders(2,i)+F(2,2)*ders(1,i)   F(3,1)*ders(2,i)+F(3,2)*ders(1,i);
                    (F(1,3)*ders(2,i) + N(i)*(k0(3)*F(2,2)-k0(2)*F(3,2)))  (F(2,3)*ders(2,i) + N(i)*(k0(1)*F(3,2)-k0(3)*F(1,2)))   (F(3,3)*ders(2,i) + N(i)*(k0(2)*F(1,2)-k0(1)*F(2,2)));
                    (F(1,3)*ders(1,i) + N(i)*(k0(3)*F(2,1)-k0(2)*F(3,1)))  (F(2,3)*ders(1,i) + N(i)*(k0(1)*F(3,1)-k0(3)*F(1,1)))   (F(3,3)*ders(1,i) + N(i)*(k0(2)*F(1,1)-k0(1)*F(2,1))) ];
    
    
                % Compute BN_6_q for all N=1:nn & q=1:6
                BN_6_q(i*3-2:i*3, :) = ders(1,i) * deps0k0_dq - N(i) * dkdq_cross_F1;
                BN_5_q(i*3-2:i*3, :) = ders(2,i) * deps0k0_dq - N(i) * dkdq_cross_F2;
                BN_3_q(i*3-2:i*3, :) = -N(i) * (dkdq_cross_F3 + cross(repmat(k0, 1, 6), deps0k0_dq));
            end
    
            % Compute R_y_partial_all
            R_y_partial_all = BN' * dtan * E_v_0_q + (BN_3_q * pk2(3) + BN_5_q * pk2(5) + BN_6_q * pk2(6));
    
            % Combine with S_reduced and integrate
            stress_resultant = stress_resultant + fac * H' * S_red;
    
            % Integrate the sensitivities residual vector
            R_y_all(sctrB, :) = R_y_all(sctrB, :) + fac * R_y_partial_all;
        end
    end
    
    % Pad R_y_all
    R_y_all_padded = zeros(ndofs+6, 6);
    R_y_all_padded(1:ndofs, :) = R_y_all;

    % Solve u_q for q=1:6
    uy_all = K \ (-R_y_all_padded);

    % Extract only sensitivity of displacement solution u
    uy_all = uy_all(1:ndofs, :);

    % Extract n0 and m0
    n0 = stress_resultant(1:3)';
    m0 = stress_resultant(4:6)';

    % Second integration to compute the C0 Entries in one integration run
    % according to (eq. 141)
    C0 = zeros(6,6); 

    for el = 1:mesh.nElems                % loop over elements
        sctr = mesh.elNodeCnt(el,:);       % element control points index
        elDoma = mesh.elDoma(el,:);        % element parametric domain
        elCpts0 = mesh.initcoords(sctr,:); % initial coordinates of el cont points
        nn = length(sctr);                % number of control points for each element
        nnElem = nn*dof;                  % total dof for each element
        sctrB = zeros(1, nnElem);

        for i = 1:dof
            sctrB(i:dof:nnElem) = dof*(sctr-1) + i;  % displacement in i-th direction
        end

        elDisp = u(sctrB);
        elDisp = reshape(elDisp, dof, nn);
        elCpts(:,1:3)=elCpts0(:,1:3)+elDisp'; %here we actualize x+du

        % Element Displacement Sensitivity over q (always last dimension)
        elY = reshape(uy_all(sctrB, :), dof, nn, 6);

        for ipt = 1:size(gp,1)            % loop over integration points
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

            dx_alpha = elDisp * ders3D';
            F = def_gradient(eps0, k0, x, dx_alpha);

            if (det(F))<0 % Check determinate 
                warning('det(F) is negative, error')
            end

            % Compute the material response
            if (mat.index >= 10 && mat.index < 20) % Formulation with PK1
                error("This solution scheme for the beam effects and stiffness is only usable with a PK2 material formulation!")
            elseif (mat.index >= 110 && mat.index < 120) % Formulation with PK2
                [pk2, dtan] = material_CSWP_PK2_hyperelasticity(dof, mat, F);
            end

            % Extract the relevant stress components for PK2_reduced 
            % (eq. 127)
            S_red = [pk2(6); pk2(5); pk2(3)];

            % Extract the relevant columns of the reduced material stiffness matrix
            % (eq. 124)
            C_red = dtan(:, [6, 5, 3]);

            % Compute the E_voigt components over q=1:6
            % (eq. 130, 123)
            H_eps0 = F';
            dkdq_cross_x = [0,x(3),-x(2); -x(3), 0, x(1); x(2), -x(1), 0];
            deps0k0_dq = [eye(3,3), dkdq_cross_x];
            H_k0 = F' * dkdq_cross_x;
            H = [H_eps0, H_k0];

            % Compute BN matrix entries (eq. 88)
            BN = zeros(6, nn*3);

            % Indices for blocks
            idx1 = 1:3:3*nn;
            idx2 = 2:3:3*nn;
            idx3 = 3:3:3*nn;

            C3 = [ (k0(3)*F(2,3)-k0(2)*F(3,3)), (k0(1)*F(3,3)-k0(3)*F(1,3)), (k0(2)*F(1,3)-k0(1)*F(2,3)) ];
            C5 = [ (k0(3)*F(2,2)-k0(2)*F(3,2)), (k0(1)*F(3,2)-k0(3)*F(1,2)), (k0(2)*F(1,2)-k0(1)*F(2,2)) ];
            C6 = [ (k0(3)*F(2,1)-k0(2)*F(3,1)), (k0(1)*F(3,1)-k0(3)*F(1,1)), (k0(2)*F(1,1)-k0(1)*F(2,1)) ];


            % Row 1
            BN(1, idx1) = F(1,1)*ders(1,:); 
            BN(1, idx2) = F(2,1)*ders(1,:); 
            BN(1, idx3) = F(3,1)*ders(1,:);

            % Row 2
            BN(2, idx1) = F(1,2)*ders(2,:); 
            BN(2, idx2) = F(2,2)*ders(2,:); 
            BN(2, idx3) = F(3,2)*ders(2,:);

            % Row 3
            BN(3, idx1) = C3(1)*N;          
            BN(3, idx2) = C3(2)*N;          
            BN(3, idx3) = C3(3)*N;

            % Row 4
            BN(4, idx1) = F(1,1)*ders(2,:) + F(1,2)*ders(1,:);
            BN(4, idx2) = F(2,1)*ders(2,:) + F(2,2)*ders(1,:);
            BN(4, idx3) = F(3,1)*ders(2,:) + F(3,2)*ders(1,:);

            % Row 5
            BN(5, idx1) = C5(1)*N + F(1,3)*ders(2,:);
            BN(5, idx2) = C5(2)*N + F(2,3)*ders(2,:);
            BN(5, idx3) = C5(3)*N + F(3,3)*ders(2,:);

            % Row 6
            BN(6, idx1) = F(1,3)*ders(1,:) + C6(1)*N;
            BN(6, idx2) = F(2,3)*ders(1,:) + C6(2)*N;
            BN(6, idx3) = F(3,3)*ders(1,:) + C6(3)*N;

            % Compute explicit strain derivatives
            dkdq_cross_x = [0,x(3),-x(2); -x(3), 0, x(1); x(2), -x(1), 0];
            deps0k0_dq = [eye(3,3), dkdq_cross_x];

            % Compute the E,q^{0} operator
            % (eq. 135)
            E_v_0_q = zeros(6,6);

            % Only fill the 3rd, 5th and 6th layers
            E_v_0_q(3, :) = F(:, 3)' * deps0k0_dq;
            E_v_0_q(5, :) = F(:, 2)' * deps0k0_dq;
            E_v_0_q(6, :) = F(:, 1)' * deps0k0_dq;

            % Compute the u_{,y}-related component of EE_v from uy_all
            % (eq. 137)
            EE_v = E_v_0_q + BN * uy_all(sctrB, :);

            % Determine local u_q and u_q_alpha
            u_q = sum(N .* elY, 2);
            u_q_alpha1 = sum(ders(1, :) .* elY, 2);
            u_q_alpha2 = sum(ders(2, :) .* elY, 2);

            % Assemble F_q following eq. 55
            k_cross_uq = cross(repmat(k0,1,1,6), u_q);
            F_q_imp = [u_q_alpha1, u_q_alpha2, k_cross_uq];
            F_q_exp = zeros(3,3,6);
            F_q_exp(:, 3, :) = deps0k0_dq;

            F_q = F_q_imp + F_q_exp;
            F_q_T = permute(F_q, [2, 1, 3]);

            % Compute a_p Term following eq. 13 in (1)
            a_p = [eye(3,3), cross(eye(3,3), repmat(x, 1, 3))];

            % Compute the second derivative of h for all q=1:6 and p=1:6
            % Add first component of HH
            HH1 = zeros(3,6,6);
            for qqi = 1:6
                HH1(:, :, qqi) = F_q_T(:,:,qqi) * a_p;
            end

            % Add second component of HH (eq. 140)
            HH2 = zeros(3,6,6);
            for ppi = 4:6 % First 3 entries in dk0dq are 0 anyways
                for qqi = 1:6
                    HH2(:, ppi, qqi) = F' * cross(dk0dq(:, ppi), u_q(:,1,qqi));
                end
            end
            
            % Combine into HH expression
            HH = HH1 + HH2;

            % Finally assemble into Beam Stiffness Components (eq. 141)
            C0_partial_geo = reshape(pagemtimes(repmat(S_red', 1,1,6), HH), 6, 6);
            C0_partial_mat = H' * C_red' * EE_v;
            C0_partial = C0_partial_mat + C0_partial_geo;
            C0 = C0 + fac * C0_partial;
        end
    end

end