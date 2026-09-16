function [omesh] = build_omesh(geo)
% This function constructs an O-mesh consisting of 4 geometrical arc-segments surrounding a square
% center segment.
% The resulting Multi-Mesh object has the same properties as a single mesh
% object with the addition of the field "submeshes", which contains the
% individual meshes (in this case 5)
% Input:
    % geo           - Employed IGA Geometry 
% Output:
    % omesh         - combined mesh object
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

n_sub = 5; % Number of Patches to actually connect into one
% May be changed to only use a subset of the patches. 
% Attention: Expect translation in the solution in these cases.


% build iga mesh structure
num_meshes = min(length(geo), n_sub);
mesh = cell(1,num_meshes);
for i = 1:num_meshes
    mesh{1,i} = build_iga_mesh( geo{1,i} );   
    mesh{1,i}.gloElNodeCnt = mesh{1,i}.elNodeCnt;
    mesh{1,i}.nodeNum = zeros(mesh{1,i}.nCptsV, mesh{1,i}.nCptsU);
end

% assuming all circle segments have the same rows and columns
n_rows = mesh{1,1}.nCptsV; % Num Points along V-Axis
n_cols = mesh{1,1}.nCptsU; % Num Points along U-Axis
p = mesh{1,1}.p;
q = mesh{1,1}.q;
pq = (p+1)*(q+1);

% Overlap of El.-Cpts for number of connected elements
rem = [0, 
    n_rows,
    2*n_rows,
    4*n_rows,
    4*n_rows + 4*(n_cols-1)];

n_points_total = num_meshes * (mesh{1,1}.nCpts) - rem(num_meshes);% Total of unique points
n_elems_total = num_meshes * mesh{1,1}.nElems;

% node numbering (top left start, row-wise)
for k = 1:num_meshes
    for j = 1:n_rows % Radius Count
        for i = 1:n_cols % Radial Count
            mesh{1,k}.nodeNum(j,i) = (j-1)*(n_cols)+i;
        end
    end
end


globCoords = zeros(n_points_total, 4);
globCoords(1:(n_rows*n_cols) , :) = mesh{1,1}.coords; % assign global coordinates of the first patch

if num_meshes >= 2
    %  node numbering & globCoords assignment for patch 2
    nodal_num = n_rows*n_cols;
    mesh{1,2}.nodeNum(1,:) = flip(mesh{1,1}.nodeNum(:,end)); % flipped last Column of Patch1 = First row of Patch2
    for j = 2:n_rows % Rows going up (+Y)
        for i = 1:n_cols % Columns going right (+X)
            bb = (j-1)*(n_rows)+i;    % Local Node Index
            nodal_num = nodal_num + 1;          % global Node Index
            
            mesh{1,2}.nodeNum(j,i) = nodal_num; % Assign correct GLOBAL Node Index
            globCoords(nodal_num, :) = mesh{1,2}.coords(bb, :); % Assign Local coords as global coords (Index-Shift)
        end
    end
end

if num_meshes >= 3
    % node Numbering & globCoords assignment for patch 3
    mesh{1,3}.nodeNum(:, end) = mesh{1,2}.nodeNum(end, :); % Last Row of Patch2 = Last Column of Patch3
    for j = 1:n_rows % Rows going up (+Y)
        for i = 1:n_cols-1 % Columns going right (+X)
            bb = (j-1)*(n_rows)+i;    % Local Node Index
            nodal_num = nodal_num + 1;          % global Node Index
            
            mesh{1,3}.nodeNum(j,i) = nodal_num; % Assign correct GLOBAL Node Index
            globCoords(nodal_num, :) = mesh{1,3}.coords(bb, :); % Assign Local coords as global coords (Index-Shift)
        end
    end
end

if num_meshes >= 4
    % node Numbering & globCoords assignment for patch 4
    mesh{1,4}.nodeNum(1, :) = mesh{1,1}.nodeNum(:, 1); % First Row of Patch4 = First Column of Patch1
    mesh{1,4}.nodeNum(end, :) = flip(mesh{1,3}.nodeNum(:, 1)); % last row of Patch4 = Last Row of Patch3
    for j = 2:n_rows-1 % Rows going up (+Y)
        for i = 1:n_cols % Columns going right (+X)
            bb = (j-1)*(n_rows)+i;    % Local Node Index
            nodal_num = nodal_num + 1;          % global Node Index
            
            mesh{1,4}.nodeNum(j,i) = nodal_num; % Assign correct GLOBAL Node Index
            globCoords(nodal_num, :) = mesh{1,4}.coords(bb, :); % Assign Local coords as global coords (Index-Shift)
        end
    end
end

if num_meshes >= 5
    % node Numbering & globCoords assignment of Patch5
    mesh{1,5}.nodeNum(1, :) = mesh{1,1}.nodeNum(end, :); % First row of Patch5 = last row of Patch1
    mesh{1,5}.nodeNum(end, :) = mesh{1,3}.nodeNum(1, :); % last row of Patch5 = first row of Patch3
    mesh{1,5}.nodeNum(:, 1) = mesh{1,4}.nodeNum(:, end); % First Column of Patch5 = Last column of Patch4
    mesh{1,5}.nodeNum(:, end) = mesh{1,2}.nodeNum(:, 1); % Last Column of Patch5 = First Column of Patch2
    
    for j = 2:n_rows-1 % Rows going up (+Y)
        for i = 2:n_cols-1 % Columns going right (+X)
            bb = (j-1)*(n_rows)+i;    % Local Node Index
            nodal_num = nodal_num + 1;          % global Node Index
            
            mesh{1,5}.nodeNum(j,i) = nodal_num; % Assign correct GLOBAL Node Index
            globCoords(nodal_num, :) = mesh{1,5}.coords(bb, :); % Assign Local coords as global coords (Index-Shift)
        end
    end
end

%% Global Connectivity

global_node_id = 0;
global_ElNodeCnt = zeros(n_elems_total, (p+1)*(q+1));
global_elDoma = zeros(n_elems_total, 4);

for kk = 1:num_meshes
    % preallocate global element connectivity matrix for the current path
    mesh{1,kk}.gloElNodeCnt = zeros(mesh{1,kk}.nElems, mesh{1,kk}.nElemCpts);
    
    el_idx = 0;
    
    % Loop over all elements in the patch
    for j = 1:mesh{1,kk}.nElemV
        for i = 1:mesh{1,kk}.nElemU
            global_node_id = global_node_id + 1; % Node index GLOBALLY
            el_idx = el_idx + 1; % Node index in THIS Mesh
            
            % Extract the neighboring nodes (based on original node
            % numbering)
            local_window = mesh{1,kk}.nodeNum(j:j+q, i:i+p);
            
            % Reduce window matrix into a row
            mesh{1,kk}.gloElNodeCnt(el_idx, :) = reshape(local_window', 1, []);
            
            global_ElNodeCnt(global_node_id, :) = reshape(local_window', 1, []);
            global_elDoma(global_node_id, :) = mesh{1, kk}.elDoma(el_idx, :);
            
        end
    end
end


%% Finalize into Mesh form
% Exemplary Format:
    %x          dim: 2
    %x            p: 3
    %x            q: 3
    %x       uKnots: [0 0 0 0 0.2500 0.5000 0.7500 1 1 1 1]
    %x       vKnots: [0 0 0 0 0.2500 0.5000 0.7500 1 1 1 1]
    %x       nCptsU: 7
    %x       nCptsV: 7
    %x        nCpts: 49
    %x       coords: [49×4 double]
    %x   initcoords: [49×4 double]
    %       nElemU: 4
    %       nElemV: 4
    %x       nElems: 16
    %x    nElemCpts: 16
    %x    elNodeCnt: [16×16 double]
    %x       elDoma: [16×4 double]
    %x gloElNodeCnt: [16×16 double]
    %      nodeNum: [7×7 double]


omesh.dim = mesh{1,1}.dim;
omesh.p = mesh{1,1}.p;
omesh.q = mesh{1,1}.q;
%omesh.k = mesh{1,1}.p; % Copy to use for 3D Mesh Visualization options
omesh.nCpts = n_points_total;
omesh.coords = globCoords;
omesh.initcoords = globCoords;
omesh.nElems = n_elems_total;
omesh.elDoma = global_elDoma;
omesh.nElemCpts = pq;
omesh.elNodeCnt = global_ElNodeCnt;
omesh.uKnots = mesh{1,1}.uKnots;
omesh.vKnots = mesh{1,1}.vKnots;
%omesh.wKnots = zeros(size(mesh{1,1}.vKnots)); % Copy to use for 3D Mesh Visualization options
omesh.nCptsU = n_cols;
omesh.nCptsV = n_rows;
omesh.submeshes = mesh;

end