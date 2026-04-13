function [ Kglob, Rglob ] = globalstiffness_CSWP_PK2( eltype, geo, mesh, mat, u , curtime,eps0,k0)
% Computes and assemlbes the global stiffness Matrix and global residuals
% vector for the cross-sectional warping problem using a hyperelastic
% material, formulated using the second Piola-Kirchhoff-Stress
% Input:
    % eltype    - (Int) Element type identifier, 30 for CSWP
    % geo       - IGA Geometry object as foound in "geometries"
    % mesh      - Mesh object, see  "build_iga_mesh(geo)"
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


if eltype == 30 %More element types cmay come in the future
    dof = 3;
end
gp_x = mesh.p+1;        % number of integration points in x-direction
gp_y = mesh.q+1;        % number of integration points in y-direction
[gp, wgt] = gauss_quadrature(gp_x, gp_y);   % calculate integration points and its weights

ndofs = dof * mesh.nCpts;      % total dofs

K = sparse(ndofs,ndofs);             % reserve stiffness matrix
C = sparse(3,ndofs);
B = sparse(3,ndofs);
R = zeros(ndofs,1);               % reserve residual matrix
Rlambda = zeros(3,1);
Rmu = zeros(3,1);
eps0 = curtime*eps0;
k0 = curtime*k0;

%This is the skew symmetric cross-product matrix
k0x = [    0   -k0(3)   k0(2);
           k0(3)     0   -k0(1);
          -k0(2)   k0(1)     0  ];  

for el = 1:mesh.nElems                % loop over elements
    sctr = mesh.elNodeCnt(el,:);       % element control points index
    elDoma = mesh.elDoma(el,:);        % element parametric domain
    %elCpts = mesh.coords(sctr,:);     % coordinates of element control points
    elCpts0 = mesh.initcoords(sctr,:); % initial coordinates of el cont points
    nn = length(sctr);                % number of control points for each element
    nnElem = nn*dof;                  % dof for each element
    sctrB = zeros(1, nnElem); 
    RE2 = zeros(nnElem,1);
    EB = zeros(3,nnElem);
    EC = zeros(3,nnElem);
    for i = 1:dof
        sctrB(i:dof:nnElem) = dof*(sctr-1) + i;  % displacement in i-th direction
    end
    
    elDisp = u(sctrB);
    elDisp = reshape(elDisp, dof, nn);
    lambda = u(end-5:end-3);
    mu = u(end-2:end);
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
        e = eye(3);
        x = N.*elCpts(:,1:dof)';
        x = sum(x,2);
        x0 = N.*elCpts0(:,1:dof)';
        x0 = sum(x0,2);
        [M,m,Theta] = GeometricalTerms(x0,x,mu,e);

        dx_alpha = elDisp * ders3D';
        F = def_gradient(eps0, k0, x, dx_alpha);

        % Retrieve PK2 material response as PK2 Stress and dtangent 
        [ stress, dtan ] = material_CSWP_PK2_hyperelasticity( dof, mat, F );
        BN = zeros(6, nn*3);
        BG = zeros(6,nn,nn);
  
        for i = 1:nn
            % Corresponds to Eq. 52 in (1)
            BN(:,i*3-2:i*3) = [ F(1,1)*ders(1,i)     F(2,1)*ders(1,i)      F(3,1)*ders(1,i);
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
            RE2(3*i-2:3*i) = N(i)*(lambda+M*mu);
            EB(:,3*i-2:3*i) = N(i)*e;
            EC(:,3*i-2:3*i) = N(i)*M';
        end
        
        G = zeros(nnElem);
        % Corresponds to Eq. (71)
        for ii=1:nn
            for jj=1:nn
                G(ii*3-2:ii*3,jj*3-2:jj*3) = eye(3)*(stress(1)*BG(1,ii,jj)+stress(2)*BG(2,ii,jj)+stress(4)*BG(4,ii,jj))...
                    +stress(3)*BG(3,ii,jj)*(k0x*k0x)+(stress(5)*BG(5,ii,jj)+stress(6)*BG(6,ii,jj))*k0x;
            end
        end
        
        RE = BN'*stress;
        R(sctrB) = R(sctrB) + fac*(RE+RE2);
        Rlambda = Rlambda + fac*x;
        Rmu = Rmu + fac*m;

        EK = BN'*dtan*BN + G + kron(N(:)*N(:)',Theta);
        
        % assemble global stiffness
        K(sctrB,sctrB) = K(sctrB,sctrB) + fac*(EK);    
        B(:,sctrB) = B(:,sctrB) + fac*EB;
        C(:,sctrB) = C(:,sctrB) + fac*EC;
    end
end
% Add additional lagrange multiplier components
Kglob = [K, B',C';B zeros(3,6);C zeros(3,6)];
Rglob = [R;Rlambda;Rmu];

end