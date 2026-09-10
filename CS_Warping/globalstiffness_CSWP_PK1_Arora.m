function [Kglob, Rglob] = globalstiffness_CSWP_PK1_Arora(eltype, geo, mesh, mat, u, curtime, eps0, k0)
% This function computes and assembles the global stiffness matrix and
% residual vector for the cross-sectional warping problem using the PK1
% formulation of Arora et al. It is the PK1 counterpart of
% globalstiffness_CSWP_PK2.
% Input:
    % eltype    - (Int) Element type identifier, 30 for CSWP
    % geo       - IGA Geometry object as found in "geometries"
    % mesh      - Mesh object, see "build_iga_mesh(geo)"
    % mat       - (Struct) Containing the material properties, see "default_mat()"
    % u         - Displacement solution vector
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
% warping problem for hyperelastic beams: A compact formulation in
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

if eltype == 30
    dof = 3;
end

gp_x = mesh.p + 1;        % number of integration points in x-direction
gp_y = mesh.q + 1;        % number of integration points in y-direction
[gp, wgt] = gauss_quadrature(gp_x, gp_y);

ndofs = dof * mesh.nCpts;

K = sparse(ndofs, ndofs);
C = sparse(3, ndofs);
B = sparse(3, ndofs);

R = zeros(ndofs, 1);
Rlambda = zeros(3, 1);
Rmu = zeros(3, 1);

eps0 = curtime * eps0;
k0   = curtime * k0;

% Skew matrix: k0x*a = k0 x a
k0x = [    0    -k0(3)   k0(2);
        k0(3)      0    -k0(1);
       -k0(2)   k0(1)      0  ];

e = eye(3);

for el = 1:mesh.nElems
    sctr = mesh.elNodeCnt(el,:);       % element control points index
    elDoma = mesh.elDoma(el,:);        % element parametric domain
    elCpts0 = mesh.initcoords(sctr,:); % initial coordinates of element control points

    nn = length(sctr);                 % number of control points for each element
    nnElem = nn * dof;                 % dofs for each element

    sctrB = zeros(1, nnElem);

    for i = 1:dof
        sctrB(i:dof:nnElem) = dof * (sctr - 1) + i;
    end

    elDisp = u(sctrB);
    elDisp = reshape(elDisp, dof, nn);

    lambda = u(end-5:end-3);
    mu     = u(end-2:end);

    elCpts(:,1:3) = elCpts0(:,1:3) + elDisp.';

    for ipt = 1:size(gp,1)
        pt = gp(ipt,:);
        wt = wgt(ipt);

        gauPts = parameter_gauss_mapping(elDoma, pt);
        j1 = jacobian_gauss_mapping(elDoma);

        [N, ders] = nurbs_derivatives(gauPts, geo, mesh);

        jmatrix = ders * elCpts0(:,1:dof-1); % mapping is 2D
        j2 = det(jmatrix);

        ders = jmatrix \ ders;

        fac = j1 * j2 * wt;

        ders3D = zeros(3, size(elCpts,1));
        ders3D(1:2,:) = ders;

        x = N .* elCpts(:,1:dof).';
        x = sum(x,2);

        x0 = N .* elCpts0(:,1:dof).';
        x0 = sum(x0,2);

        [M, m, Theta] = GeometricalTerms(x0, x, mu, e);

        dx_alpha = elDisp * ders3D.';
        F = def_gradient(eps0, k0, x, dx_alpha);

        % PK1 material response
        [pk1, A] = material_CSWP_hyperelasticity(dof, mat, F);

        % Precompute tangent slices at the current Gauss point
        A_i3_j3 = squeeze(A(:,3,:,3));

        A_i3_j1 = squeeze(A(:,3,:,1));
        A_i3_j2 = squeeze(A(:,3,:,2));

        A_i1_j3 = squeeze(A(:,1,:,3));
        A_i2_j3 = squeeze(A(:,2,:,3));

        A_i1_j1 = squeeze(A(:,1,:,1));
        A_i1_j2 = squeeze(A(:,1,:,2));
        A_i2_j1 = squeeze(A(:,2,:,1));
        A_i2_j2 = squeeze(A(:,2,:,2));

        % PK1 stress vector, column-wise
        Pvec = pk1(:);

        RE  = zeros(nnElem, 1);
        RE2 = zeros(nnElem, 1);
        EK  = zeros(nnElem, nnElem);

        EB = zeros(3, nnElem);
        EC = zeros(3, nnElem);

        % Residual: Arora Eq. (39)
        % H_J = N_J (P e3 x k0) + P e_alpha N_{J,alpha}
        % Since P e3 x k0 = pk1(:,3) x k0 = -k0x * pk1(:,3)
        for J = 1:nn
            rows = 3*(J-1) + (1:3);

            NJ  = N(J);
            NJ1 = ders(1,J);
            NJ2 = ders(2,J);

            QJ = [ ...
                NJ1 * e, ...
                NJ2 * e, ...
               -NJ  * k0x ];

            RE(rows) = QJ * Pvec;

            % Constraint residual and Lagrange multiplier couplings
            RE2(rows) = NJ * (lambda + M * mu);

            EB(:,rows) = NJ * e;
            EC(:,rows) = NJ * M.';
        end

        % Stiffness: optimized Arora A_JI block assembly, Eq. (37)
        for J = 1:nn
            rows = 3*(J-1) + (1:3);

            NJ  = N(J);
            NJ1 = ders(1,J);
            NJ2 = ders(2,J);

            for I = 1:nn
                cols = 3*(I-1) + (1:3);

                NI  = N(I);
                NI1 = ders(1,I);
                NI2 = ders(2,I);

                % Combined terms T1 + T2
                T12 = -k0x * ( ...
                          A_i3_j3 * k0x * NI ...
                        + A_i3_j1 * NI1 ...
                        + A_i3_j2 * NI2 ) * NJ;

                % Combined terms T3 + T4
                T34 = NJ1 * ( ...
                          A_i1_j3 * k0x * NI ...
                        + A_i1_j1 * NI1 ...
                        + A_i1_j2 * NI2 ) ...
                    + NJ2 * ( ...
                          A_i2_j3 * k0x * NI ...
                        + A_i2_j1 * NI1 ...
                        + A_i2_j2 * NI2 );

                EK(rows,cols) = T12 + T34;
            end
        end

        % Constraint stiffness: block (J,I) = N(J) * N(I) * Theta
        EK = EK + kron(N(:) * N(:).', Theta);

        % Assemble residuals
        R(sctrB) = R(sctrB) + fac * (RE + RE2);
        Rlambda = Rlambda + fac * x;
        Rmu     = Rmu     + fac * m;

        % Assemble stiffness and Lagrange multiplier coupling blocks
        K(sctrB,sctrB) = K(sctrB,sctrB) + fac * EK;
        B(:,sctrB) = B(:,sctrB) + fac * EB;
        C(:,sctrB) = C(:,sctrB) + fac * EC;
    end
end

% Add Lagrange multiplier components
Kglob = [K, B', C';
         B, zeros(3,6);
         C, zeros(3,6)];

Rglob = [R; Rlambda; Rmu];

end