function [ Kglob, Rglob ] = globalstiffness_CSWP_PK1_Arora( eltype, geo, mesh, mat, u, curtime, eps0, k0 )
% Computes and assembles the global stiffness matrix and residual vector
% for the cross-sectional warping problem using the PK1 formulation of
% Arora et al. implemented in an optimized blockwise manner.
%
% This function is the PK1 counterpart of globalstiffness_CSWP_PK2.
%
% Required material routine:
%   [ pk1, A ] = material_CSWP_hyperelasticity( dim, mat, F );
%
% where:
%   pk1          : 3 x 3 first Piola-Kirchhoff stress
%   A(i,A,j,B)  : dP_iA / dF_jB, size 3 x 3 x 3 x 3


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

        % ------------------------------------------------------------
        % PK1 material response
        % pk1 : 3 x 3
        % A   : 3 x 3 x 3 x 3, A(i,A,j,B) = dP_iA/dF_jB
        % ------------------------------------------------------------
        [pk1, A] = material_CSWP_hyperelasticity(dof, mat, F);

        % ------------------------------------------------------------
        % Precompute tangent slices once per Gauss point
        % ------------------------------------------------------------
        A_i3_j3 = squeeze(A(:,3,:,3));

        A_i3_j1 = squeeze(A(:,3,:,1));
        A_i3_j2 = squeeze(A(:,3,:,2));

        A_i1_j3 = squeeze(A(:,1,:,3));
        A_i2_j3 = squeeze(A(:,2,:,3));

        A_i1_j1 = squeeze(A(:,1,:,1));
        A_i1_j2 = squeeze(A(:,1,:,2));
        A_i2_j1 = squeeze(A(:,2,:,1));
        A_i2_j2 = squeeze(A(:,2,:,2));

        % PK1 stress vector, column-wise:
        % Pvec = [pk1(:,1); pk1(:,2); pk1(:,3)]
        Pvec = pk1(:);

        RE  = zeros(nnElem, 1);
        RE2 = zeros(nnElem, 1);
        EK  = zeros(nnElem, nnElem);

        EB = zeros(3, nnElem);
        EC = zeros(3, nnElem);

        % ------------------------------------------------------------
        % Residual: Arora Eq. (39)
        %
        % H_J = N_J (P e3 x k0) + P e_alpha N_{J,alpha}
        %
        % Since P e3 x k0 = pk1(:,3) x k0 = -k0x * pk1(:,3)
        % ------------------------------------------------------------
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

        % ------------------------------------------------------------
        % Stiffness: optimized Arora A_JI block assembly, Eq. (37)
        % ------------------------------------------------------------
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

        % ----------------------------------------------------
        % Combined terms T1 + T2:
        %
        % T1 = -k0x * A_i3_j3 * k0x * NI * NJ
        % T2 = -k0x * (A_i3_j1*NI1 + A_i3_j2*NI2) * NJ
        %
        % Common factor: -k0x * (...) * NJ
        % ----------------------------------------------------
        T12 = -k0x * ( ...
                  A_i3_j3 * k0x * NI ...
                + A_i3_j1 * NI1 ...
                + A_i3_j2 * NI2 ) * NJ;

        % ----------------------------------------------------
        % Combined terms T3 + T4:
        %
        % T3 = (NJ1*A_i1_j3 + NJ2*A_i2_j3) * k0x * NI
        %
        % T4 = NJ1*(A_i1_j1*NI1 + A_i1_j2*NI2)
        %    + NJ2*(A_i2_j1*NI1 + A_i2_j2*NI2)
        %
        % Common factors: NJ1 and NJ2
        % ----------------------------------------------------
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

        % Constraint stiffness:
        % block (J,I) = N(J) * N(I) * Theta
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