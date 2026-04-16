function [n0_entry] = beam_forces_entry(geo, mesh, mat, eps0, k0, u, ip)
% Computate a beam force entry (as outlined in "Numerische Methoden zur 
% Modellierung elastoplastischer Balken und ihre Anwendung auf 
% periodische Gitterstrukturen", PhD Thesis by L. Herrnböck)
% Input:
	% geo   - Employed IGA Geometry 
	% mesh  - Employed mesh 
	% mat   - struct containing material parameters
    % eps0  - (3,1) vector containing the strain prescriptors
    % k0    - (3,1) vector containing the twist prescriptors
    % u     - displacement solution vector
    % ip    - selector index (corresponds to n0 in range 1:3 and m0 in range 4:6)
% Output:
    % n0_entry  - (3,1) force/moment related to selector index (in kN or kNmm)
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
% notation", DOI: 10.48550/arXiv.2604.12886 
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
[gp, wgt] = gauss_quadrature(gp_x, gp_y);   % calculate integration points and its weights

dof = 3;
count = 0;

n0_entry = 0; % Reserve beam force entry value

% Compute derivatives of v0 and k0 for p and q -> selection vectors
[deps0_dp, dk0_dp] = derivatives_vk(eps0, k0, ip);

% Prenotate the unit tensor
e = eye(3);

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

        dx_alpha = elDisp * ders3D';
        F = def_gradient(eps0, k0, x, dx_alpha);

        if (det(F))<0 % Check determinate 
            warning('det(F) is negative, error')
        end
        
        if (mat.index >= 10 && mat.index < 20) % Formulation with PK1
            [ pk1, ~ ] = material_CSWP_hyperelasticity(dof, mat, F);
        elseif (mat.index >= 110 && mat.index < 120) % Formulation with PK2
            [pk2, ~] = material_CSWP_PK2_hyperelasticity(dof, mat, F);

            % Compute PK1 from PK2 by transforming PK2 back to matrix-form
            pk1 = F * to_tensor(pk2, "stress");
        end

        % Compute the beam force entry
        term1 = deps0_dp * e(:,3)' + cross(dk0_dp, x) * e(:,3)';
        n0_entry_partial = tensorprod(pk1, term1, "all");
            
        % Update n0_entry at each gauss point
        n0_entry = n0_entry + fac * n0_entry_partial;
    end
end