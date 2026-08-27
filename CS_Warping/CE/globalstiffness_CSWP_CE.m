function [ Kglob, Rglob ] = globalstiffness_CSWP_CE( eltype, geo, mesh, mat, v, c0, curtime,eps0,k0, dt)
% Computes and assemlbes the global stiffness Matrix and global residuals
% vector for the cross-sectional warping problem using a hyperelastic
% material, formulated using the second Piola-Kirchhoff-Stress
% Input:
    % eltype    - (Int) Element type identifier, 30 for CSWP
    % geo       - IGA Geometry object as foound in "geometries"
    % mesh      - Mesh object, see  "build_iga_mesh(geo)"
    % mat       - (Struct) Containing the material properties, see "default_mat()"
    % v         - Chemo-displacement solution vector
    % curtime   - Current time step
    % eps0      - (3,1) vector containing the strain prescriptors
    % k0        - (3,1) vector containing the twist prescriptors
% Output:
    % Kglob     - (n,n) Matrix of global stiffness entries
    % Rglob     - (n,1) Vector of global residual entries
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
% (3) A. Shafqat, O. Weeger, B. Xu, "A robust finite strain isogeometric 
% solid-beam element", Computer Methods in Applied Mechanics and 
% Engineering, DOI: 10.48550/arXiv.2312.07124
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


if eltype == 30 %More element types cmay come in the future
    dof = 3;
end
gp_x = mesh.p+1;        % number of integration points in x-direction
gp_y = mesh.q+1;        % number of integration points in y-direction
[gp, wgt] = gauss_quadrature(gp_x, gp_y);   % calculate integration points and their weights

npts = mesh.nCpts;              % total control points
ndofs = dof * mesh.nCpts;       % total dofs (3 * npts)


% Allocate space for the solution dimensions u, c, mu
K_uu = sparse(ndofs,ndofs);             % reserve stiffness matrices
K_uc = sparse(ndofs, npts);

K_cc = zeros(npts, npts);
K_cmu = zeros(npts, npts);

K_muu = zeros(npts, ndofs);
K_muc = zeros(npts, npts);
K_mumu = zeros(npts, npts);


R_u = zeros(ndofs,1);               % reserve residual matrices
R_c = zeros(npts, 1);
R_mu = zeros(npts, 1);


% Allocate space for the lagrange multipliers
%Rlambda = zeros(3,1);
%Rmu = zeros(3,1);
%C = sparse(3,ndofs);
%B = sparse(3,ndofs);

R_lam1 = zeros(3,1);
R_lam2 = zeros(3,1);
K_ulam1 = sparse(3, ndofs);
K_ulam2 = sparse(3, ndofs);

eps0 = curtime*eps0;
k0 = curtime*k0;

e = eye(3);
%This is the skew symmetric cross-product matrix
k0x = [    0   -k0(3)   k0(2);
           k0(3)     0   -k0(1);
          -k0(2)   k0(1)     0  ];  

% Disassemble solution vector v
u = v(1:ndofs);
c = v(ndofs+1:ndofs+npts);
mu = v(ndofs+npts+1:end-6);
lam1 = v(end-5:end-3);
lam2 = v(end-2:end);

% Disasselbme last converged solution vector
%u0 = v0(1:ndofs);
%c0 = v0(ndofs+1:ndofs+ncpts);
%mu0 = v0(ndofs+ncpts+1:end-6);
%lam10 = v0(end-5:end-3);
%lam20 = v0(end-2:end);

for el = 1:mesh.nElems                % loop over elements
    sctr = mesh.elNodeCnt(el,:);       % element control points index
    elDoma = mesh.elDoma(el,:);        % element parametric domain
    %elCpts = mesh.coords(sctr,:);     % coordinates of element control points
    elCpts0 = mesh.initcoords(sctr,:); % initial coordinates of el cont points
    nn = length(sctr);                % number of control points for each element
    nnElem = nn*dof;                  % dof for each element
    
    sctrB = zeros(1, nnElem); 
    
    % allocate space for element stiffness components
    

    % Allocate space for the element residual components
    %RE2 = zeros(nnElem,1);
    %EB = zeros(3,nnElem);
    %EC = zeros(3,nnElem);
    
    ER_u2 = zeros(nnElem, 1);
    ER_u3 = zeros(nnElem, 1);

    % Allocate space for the Lagrange multiplier components
    E_lam1 = zeros(3, nnElem);
    E_lam2 = zeros(3, nnElem);

    for i = 1:dof
        sctrB(i:dof:nnElem) = dof*(sctr-1) + i;  % displacement in i-th direction
    end
    
    % Element Displacement, concentration and potential
    elDisp = u(sctrB);
    elDisp = reshape(elDisp, dof, nn);
    
    elc = c(sctr);
    elc0 = c0(sctr);
    elmu = mu(sctr);
    
    %lam1 = u(end-5:end-3);
    %mu = u(end-2:end);
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
        
        % Compute Integartion point parameters x, c, mu
        x = N.*elCpts(:,1:dof)';
        x = sum(x,2);

        cipt = sum(elc' .* N, 2);
        c0ipt = sum(elc0' .* N, 2);
        muipt = sum(elmu' .* N, 2);

        x0 = N.*elCpts0(:,1:dof)';
        x0 = sum(x0,2);
        [M,m,Theta] = GeometricalTerms(x0,x,lam2,e);

        dx_alpha = elDisp * ders3D';
        % F = def_gradient(eps0, k0, x, dx_alpha);
        %F = def_gradient(eps0, k0, x, dx_alpha);  %TODO: Check if this F is applicable 
           
        [F, F_c, J_c, J_e, I1] = chemo_kin(eps0, k0, x, dx_alpha, cipt, mat.Omega);

        % Compute kinematic parameters
        %J_e = det(F);
        %C = F'*F;
        %I1 = (C(1,1) + C(2,2) + C(3,3));

        % Retrieve PK2 material response as PK2 Stress and tangent 
        %[ stress, dtan ] = material_CSWP_PK2_hyperelasticity( dof, mat, F );
        [ stress, dtan, dtan_c, dmu_dc, dmu_dC] = material_CSWP_CE_hyperelasticity(dof, mat, F, cipt, J_c, J_e, I1);
        
        % Chemical elastic potential due to elastic deformation
        mu_e_ipt = (mat.Omega * mat.elastic_modulus1)/2 * (1 - J_e^2) + ...
            (mat.Omega * mat.elastic_modulus2)/2 * (I1 - 3);


        %BN_C and BN_mu
        BN_c = [ders(1, :); ders(2, :); zeros(1, 16)];
        BN_mu = [ders(1, :); ders(2, :); zeros(1, 16)];

        % Gradient c and mu
        grad_c = sum(BN_c * c(sctr), 2);
        grad_mu = sum(BN_mu * mu(sctr), 2);

        BN = zeros(6, nn*3);
        BG = zeros(6,nn,nn);

        % Create Repeated shape fucntion vector
        N_rep3 = repelem(N,3);

        for i = 1:nn
            % Corresponds to Eq. 52 in (1)
            BN(:,i*3-2:i*3) = [ 
                F(1,1)*ders(1,i)     F(2,1)*ders(1,i)      F(3,1)*ders(1,i);
                F(1,2)*ders(2,i)     F(2,2)*ders(2,i)      F(3,2)*ders(2,i);
                N(i)*(k0(3)*F(2,3)-k0(2)*F(3,3))    N(i)*(k0(1)*F(3,3)-k0(3)*F(1,3))       N(i)*(k0(2)*F(1,3)-k0(1)*F(2,3)) ;
                F(1,1)*ders(2,i)+ F(1,2)*ders(1,i)  F(2,1)*ders(2,i)+F(2,2)*ders(1,i)   F(3,1)*ders(2,i)+F(3,2)*ders(1,i);
                (F(1,3)*ders(2,i) + N(i)*(k0(3)*F(2,2)-k0(2)*F(3,2)))  (F(2,3)*ders(2,i) + N(i)*(k0(1)*F(3,2)-k0(3)*F(1,2)))   (F(3,3)*ders(2,i) + N(i)*(k0(2)*F(1,2)-k0(1)*F(2,2)));
                (F(1,3)*ders(1,i) + N(i)*(k0(3)*F(2,1)-k0(2)*F(3,1)))  (F(2,3)*ders(1,i) + N(i)*(k0(1)*F(3,1)-k0(3)*F(1,1)))   (F(3,3)*ders(1,i) + N(i)*(k0(2)*F(1,1)-k0(1)*F(2,1))) ];

            for j=1:nn
                % Corresponds to Eq. 70 in (1)
                % Attention: kappa-related terms are reintroduced later in term G
                BG(:,i,j) = [ders(1,i)*ders(1,j);
                    ders(2,i)*ders(2,j);
                    -N(i)*N(j);                                
                    ders(1,i)*ders(2,j)+ders(2,i)*ders(1,j);
                    ders(2,i)*N(j)-N(i)*ders(2,j);            
                    ders(1,i)*N(j)-N(i)*ders(1,j)   
                    ];
            end
        end
        
        % Allocate Space for Element Stiffness Entries
        G = zeros(nnElem);
        %EK_uc = zeros(nn*3, nn);
        %EK_cc1 = zeros(nn, nn);
        %EK_cc2 = zeros(nn, nn);
        %EK_cmu = zeros(nn, nn);

        % Allocate Space for the Element Residual Entries
        %ER_c1 = zeros(nn,1);
        %ER_c2 = zeros(nn,1);
        %ER_c3 = zeros(nn,1);
        %ER_mu = zeros(nn,1);

        % Fill in the element stiffness and residual entries
        % Corresponds to Eq. (71)
        for ii=1:nn
            % Compute residual entries dependent upon only NI
            %ER_c1(ii) = N(ii) * (cipt - c0ipt) / dt;
            %ER_c2(ii) = (BN_c(:, ii)' * grad_mu) * cipt * (1-cipt);
            %ER_c3(ii) = 0; % TODO: how to implement boundary conditions for flux? 
            %ER_mu(ii) = N(ii) * (muipt - log((cipt / (1 - cipt)) - mu_e_ipt));
            %ER_u2(3*ii-2:3*ii) = N(ii)*(lam1 + M * lam2);

            %E_lam1(:,3*i-2:3*i) = N(ii)*e;
            %E_lam2(:,3*i-2:3*i) = N(ii)*M';

            % Compute Half vectorized Stiffness entries
            %EK_uc(:, ii) = BN' * dtan_c .* N(ii); % TODO: CHeck if moving from ii to jj here makes sense -> Because computation
            %EK_muu(ii, :) = -N(ii) * (dmu_dC * BN); % TODO: Check dimensions
           

            for jj=1:nn
                G(ii*3-2:ii*3,jj*3-2:jj*3) = eye(3)*(stress(1)*BG(1,ii,jj)+stress(2)*BG(2,ii,jj)+stress(4)*BG(4,ii,jj))...
                    +stress(3)*BG(3,ii,jj)*(k0x*k0x)+(stress(5)*BG(5,ii,jj)+stress(6)*BG(6,ii,jj))*k0x;

                % Compute unvectorized stiffness entries dependend on NI and NJ
                %EK_cc1(ii, jj) = (1/dt) * N(ii) * N(jj);
                %EK_cc2(:, jj) = ((1-2*cipt) * BN_c' * grad_c) * N(jj);
                %EK_cmu(ii, jj) = cipt * (1-cipt) * BN_c(:, ii)' * BN_mu(:, jj);
 
            end
        end
        
        % Compute complete vectorized Stiffness entries for u, c, mu
        dNN = (N' * N);
        cc1 = cipt * (1 - cipt);
        EK_uu = BN' * dtan * BN + G + kron(N(:)*N(:)',Theta);
        EK_uc = BN' * dtan_c .* N_rep3; % TODO: Check if this is equivalent and correct
        
        EK_cc = 1/dt * dyad(N, N) + (1 - 2*cipt) * dyad(BN_c', N); % TODO: Check if this is equivalent and correct
        EK_cmu = cc1 * (BN_c * BN_mu);% TODO: Check if this is equivalent and correct
        
        EK_muu = -dyad(N, dmu_dC * BN);% TODO: Check if this is equivalent and correct
        EK_mumu = dNN;
        EK_muc = -(1/(cc1) - dmu_dc) * dNN;

        % Compute Residual entries for u, c, mu
        ER_u2 = N*(lam1 + M * lam2);
        ER_u3 = zeros(nn, 1); % TODO: How to implement boundary traction / Force conditions?
        ER_u = BN'*stress + ER_u2 + ER_u3;
        
        ER_c1 = N * (cipt - c0ipt) / dt;
        ER_c2 = (BN_c' * grad_mu) * cc1;
        ER_c3 = zeros(nn,1); % TODO: How to implement boundary condition for concentration?
        ER_c = ER_c1 + ER_c2 + ER_c3;

        ER_mu = N * (muipt - log((cipt / (1 - cipt)) - mu_e_ipt));

        % Compute change to lagrande factors
        E_lam1 = N*e;
        E_lam2 = N*M';

        % assemble global stiffnesses
        K_uu(sctrB,sctrB) = K_uu(sctrB,sctrB) + fac * (EK_uu);
        K_uc(sctrB, sctr) = K_uc(sctrB, sctr) + fac * (EK_uc);
        %K_umu(sctrB, sctr) = K_umu(sctrB, sctr) + fac * (EK_umu);
        
        K_cc(sctr, sctr) = K_cc(sctr, sctr) + fac * (EK_cc);
        K_cmu(sctr, sctr) = K_cmu(sctr, sctr) + fac * (EK_cmu);

        K_muu(sctr, sctrB) = K_muu(sctr, sctrB) + fac * (EK_muu);
        K_muc(sctr, sctr) = Kmuc(sctr, sctr) + fac * (EK_muc);
        K_mumu(sctr, sctr) = K_mumu(sctr, sctr) + fac * (EK_mumu);

        % Assemble global residuals
        R_u(sctrB) = R_u(sctrB) + fac * (ER_u);
        R_c(sctr) = R_c(sctr) + fac * (ER_c);
        R_mu(sctr) = R_mu(sctr) + fac * (ER_mu);


        % Handle the Lagrange Multipliers
        %Rlambda = Rlambda + fac*x;
        %Rmu = Rmu + fac*m;
        %B(:,sctrB) = B(:,sctrB) + fac*EB;
        %C(:,sctrB) = C(:,sctrB) + fac*EC;
        R_lam1 = R_lam1 + fac * x;
        R_lam2 = R_lam2 + fac * m;
        K_ulam1(:, sctrB) = K_ulam1(:, sctrB) + fac * E_lam1;
        K_ulam2(:, sctrB) = K_ulam2(:, sctrB) + fac * E_lam2;
    end
end
% Add additional lagrange multiplier components

%Kglob = [K, B',C';B zeros(3,6);C zeros(3,6)];
%Rglob = [R;Rlambda;Rmu];

% B = K_ulam1
% C = K_ulam2
O1 = zeros(npts, 3);
O2 = zeros(ndofs, npts);
O3 = zeros(3,3);
Kglob = [
    K_uu,       K_uc,       O2,         K_ulam1',   K_ulam2';
    O2',        K_cc,       K_cmu,      O1,         O1;
    K_muu,      K_muc,      K_mumu,     O1,         O1;
    K_ulam1,    O1',        O1',        O3,         O3;
    K_ulam2,    O1',        O1',        O3,          O3
    ];

% Rlambda = R_lam1
% Rmu = R_lam2
Rglob = [R_u; R_c; R_mu; R_lam1; R_lam2];


end