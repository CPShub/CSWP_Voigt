function vmesh = output_visual_mesh_CSWP( fout, mat, geo, mesh, u, step, currentime, eps0, k0)
% This function generates the virtual mesh (vmesh) struct and writes its
% contents to a .msh file. 
% This variant is meant to work with results generated solving the CSWP
% Input:
    % fout          - file handle of visualized mesh
	% mat           - (Struct) containing material parameters
	% geo           - Employed IGA Geometry 
	% mesh          - Employed mesh 
    % u             - Displacement solution vector
    % step          - Current simulation step
    % currentime    - Current simulation timestamp   
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

% Construct polygon and kntcrv in Param space
num1 = 10;        
num2 = 10;
polygon = build_visual_mesh_suf( num1, num2);        % build visualized mesh
kntcrv = build_visual_knotcurve_suf( mesh.uKnots, mesh.vKnots, num1+1 ); % build visualized knot curves

% Pre-Assign vmesh
numpts = (num1+1)*(num2+1) + size(kntcrv.linpts,1);
vmesh.nodalpts = zeros(num_meshes * numpts,3);
vmesh.displacement = zeros(num_meshes* numpts,3);
vmesh.stress = zeros(num_meshes * numpts,6);
vmesh.strain = zeros(num_meshes * numpts,6);

% Iterate over all meshes
j = 0;
for m = 1:num_meshes
    sub_mesh = mesh_cell{1, m};
    %sub_geo = geo{1,5};

    elem_index = find_point_span( sub_mesh, polygon.tripts );
    
    % Iterate over all Polygon-Points
    for i=1:(num1+1)*(num2+1)
        xi = polygon.tripts(i,1);  % u sub_mesh point
        eta = polygon.tripts(i,2);  % v sub_mesh point
        
        e = elem_index(i);   % element number
        sctr = sub_mesh.elNodeCnt(e,:);     % element control points index
        exyz = sub_mesh.coords(sctr,:);  % element control points' coordinates
        nn = length(sctr);   % number of control points in the element
        nn3 = nn*3;          % degree of freedom of control points
        %nn3 = nn*2;
        
        
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
    

        % TODO: Remove after testing
        [N,ders] = nurbs_derivatives( [xi, eta],geo, mesh );
        %[N_sub, ders_sub] = nurbs_derivatives( [xi, eta], sub_geo, sub_mesh);
        
        %ders = ders_sub;
        %N = N_sub;

        jmatrix = ders*exyz(:,1:2); %Because the mapping is in 2D
        ders =  jmatrix \ ders;      
        ders3D = zeros(3,size(exyz,1));
        ders3D(1:2,:) = ders;
        x = N.*exyz(:,1:3)';
        x = sum(x,2);

        dx_alpha = edsp * ders3D';
        F = def_gradient(eps0, k0, x, dx_alpha);

        % Retrieve PK2 material response as PK2 Stress and dtangent 
        [ pk2, ~] = material_CSWP_PK2_hyperelasticity( 3, mat, F );
        ccy = pk2cauchy(pk2, F);
        
        
        % global vmesh Index j
        j = j + 1;
        vmesh.displacement(j,:) = (edsp * N')';
        vmesh.nodalpts(j,:) = N*exyz(:,1:3) + vmesh.displacement(j,1:3);   
        vmesh.stress(j,:) = ccy';
        strain = (F'*F-eye(3))/2;
        vmesh.strain(j,:) = voigt(strain)';
    end
end

for m = 1:num_meshes
    sub_mesh = mesh_cell{1, m};

    %count = (num1+1)*(num2+1);
    line_index = find_point_span( sub_mesh, kntcrv.linpts );
    kntcrv.linsub_mesh = kntcrv.linmesh + (num1+1)*(num2+1);

    % Iterate over all Knot-Curve Points
    for i = 1:size(kntcrv.linpts,1)
        %count = count+1;
        xi = kntcrv.linpts(i,1);
        eta = kntcrv.linpts(i,2);
        e = line_index(i);   % element number
        sctr = sub_mesh.elNodeCnt(e,:);     % element control points index
        exyz = sub_mesh.coords(sctr,:);  % element control points' coordinates
        nn = length(sctr);   % number of control points in the element
        nn3 = nn*3;          % degree of freedom of control points
        
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
    
        [N,ders] = nurbs_derivatives( [xi, eta],geo, mesh );
        jmatrix = ders*exyz(:,1:2); %Because the mapping is in 2D
        ders =  jmatrix \ ders;      
        ders3D = zeros(3,size(exyz,1));
        ders3D(1:2,:) = ders;
        x = N.*exyz(:,1:3)';
        x = sum(x,2);

        dx_alpha = edsp * ders3D';
        F = def_gradient(eps0, k0, x, dx_alpha);

        % Retrieve PK2 material response as PK2 Stress and dtangent 
        [ pk2, ~] = material_CSWP_PK2_hyperelasticity( 3, mat, F );
        cauchy = pk2cauchy(pk2, F);

        % Assign global vmesh Index j
        j = j + 1;
        vmesh.displacement(j,:) = edsp * N';
        vmesh.nodalpts(j,:) = N*exyz(:,1:3) + vmesh.displacement(j,1:3);   
        vmesh.stress(j,:) = cauchy';
        strain = (F'*F-eye(3))/2;
        vmesh.strain(i,:) = voigt(strain)';
    end
end

% Write to the output file
fprintf(fout,'STEP = %d, TIME = %e\n', step, currentime);
for i = 1:num_meshes*numpts
    fprintf(fout,'v %f %f %f\n', vmesh.nodalpts(i,:));
    fprintf(fout,'d %f %f %f\n', vmesh.displacement(i,:));
    fprintf(fout,'s %f %f %f %f %f %f\n', vmesh.stress(i,:));
    fprintf(fout,'t %f %f %f %f %f %f\n', vmesh.strain(i,:));
end
for m = 1:num_meshes
    last_node_id = (m-1)*size(polygon.tripts, 1);
    for i = 1:size(polygon.trimesh,1)
        fprintf(fout,'f %d %d %d\n', polygon.trimesh(i,:) + last_node_id);
    end
end
for i = 1:size(kntcrv.linmesh,1)
    fprintf(fout, 'l ');
    for j = 1:size(kntcrv.linmesh,2)
        fprintf(fout,'%d ', kntcrv.linmesh(i,j));
    end
    fprintf(fout, '\n');
end


end

