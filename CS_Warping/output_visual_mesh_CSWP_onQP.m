function vmesh = output_visual_mesh_CSWP_onQP( fout, mat, geo, mesh, u, step, currentime, eps0, k0)
% This function generates the virtual mesh (vmesh) struct by evaluating the 
% CSWP solution u at the integration points. The resolved solution is 
% written to a .msh file. 
% Input:
    % fout          - file handle of visualized mesh
	% mat           - (Struct) containing material parameters
	% geo           - Employed IGA Geometry 
	% mesh          - Employed mesh 
    % u             - Displacement solution vector
    % step          - Current simulation step
    % currentime    - Current simulation timestamp   
    % eps0          - Vector containing the strain prescriptors
    % k0            - Vector containing the twist prescriptors
% Output:
	% vmesh         - output visualized mesh structure, see "read_visual_mesh()"
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

% Check if mesh contains sub_meshes
if isfield(mesh, "submeshes")
    mesh_cell = mesh.submeshes;
    num_meshes = size(mesh.submeshes, 2);
    use_gloElNodeCnt = 1;
else
    mesh_cell = cell(1,1);
    mesh_cell{1,1} = mesh;
    num_meshes = 1;
    use_gloElNodeCnt = 0;
end

% Total number of points = Total number of gauss Quadrature Points
dof = 3;
gp_x = mesh.p+1;        % number of integration points in x-direction
gp_y = mesh.q+1;        % number of integration points in y-direction
[gp, ~] = gauss_quadrature(gp_x, gp_y);   % calculate integration points and their weights



% Pre-Assign vmesh
numpts = mesh.nElems * size(gp, 1);
vmesh.nodalpts = zeros(numpts,3);
vmesh.displacement = zeros(numpts,3);
vmesh.stress = zeros(numpts,6);
vmesh.strain = zeros(numpts,6);

% Iterate over all meshes
j = 0;
for m = 1:num_meshes
    sub_mesh = mesh_cell{1, m};
    %sub_geo = geo{1,5};

    for e = 1:sub_mesh.nElems
        sctr = sub_mesh.elNodeCnt(e,:);     % element control points index
        %exyz = sub_mesh.coords(sctr,:);  % element control points' coordinates
        nn = length(sctr);   % number of control points in the element
        nn3 = nn*3;          % degree of freedom of control points
        %nn3 = nn*2;
        elDoma = sub_mesh.elDoma(e, :);
        
        
        % Check if global globElNodeCnt should be used
        if use_gloElNodeCnt
            sctr = sub_mesh.gloElNodeCnt(e, :);
        end
        sctrB = zeros(1, nn3);      
        sctrB(1:3:nn3) = 3*sctr - 2;% displacement in x direction
        sctrB(2:3:nn3) = 3*sctr-1;  % displacement in y direction
        sctrB(3:3:end) = 3*sctr;    % displacement in z direction
        edsp = u(sctrB);
        edsp = reshape(edsp, 3, nn);

        elCpts0 = mesh.initcoords(sctr,:); % initial coordinates of el cont points
        elCpts(:,1:3)=elCpts0(:,1:3)+edsp';

        for ipt = 1:size(gp,1)
            pt = gp(ipt,:);      % reference parametric coordinates for each integration point
            gauPts = parameter_gauss_mapping( elDoma, pt );   % gauss integration mapping  
            [N,ders] = nurbs_derivatives( gauPts,geo, mesh );
            jmatrix = ders*elCpts0(:,1:dof-1); %Because the mapping is in 2D
            ders =  jmatrix \ ders;    
            ders3D = zeros(3,size(elCpts,1));
            ders3D(1:2,:) = ders;
            x = N.*elCpts(:,1:dof)';
            x = sum(x,2);
    
            dx_alpha = edsp * ders3D';
            F = def_gradient(eps0, k0, x, dx_alpha);
    

            if (mat.index >= 10 && mat.index < 20)
                % Retrieve PK1 material response as PK2 Stress and
                % transform into Cauchy Stress
                mat2 = mat;
                mat2.index = mat.index + 100;
                [ stress, ~ ] = material_CSWP_PK2_hyperelasticity( dof, mat2, F );
            elseif (mat.index >= 110 && mat.index < 120)
                % Retrieve PK2 material response as PK2 Stress and
                % transform into Cauchy Stress
                [ stress, ~ ] = material_CSWP_PK2_hyperelasticity( dof, mat, F );
            end
            ccy = pk2cauchy(stress, F);
        
            % global vmesh Index j
            j = j + 1;
            vmesh.displacement(j,:) = (edsp * N')';
            vmesh.nodalpts(j,:) = x;%N*exyz(:,1:3) + vmesh.displacement(j,1:3);   
            vmesh.stress(j,:) = ccy';
            strain = (F'*F-eye(3))/2;
            vmesh.strain(j,:) = voigt(strain)';
        end
    end
end

% Write to the output file
fprintf(fout,'STEP = %d, TIME = %e\n', step, currentime);
for i = 1:j
    fprintf(fout,'v %f %f %f\n', vmesh.nodalpts(i,:));
    fprintf(fout,'d %f %f %f\n', vmesh.displacement(i,:));
    fprintf(fout,'s %f %f %f %f %f %f\n', vmesh.stress(i,:));
    fprintf(fout,'t %f %f %f %f %f %f\n', vmesh.strain(i,:));
end
end

