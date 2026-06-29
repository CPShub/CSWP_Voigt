function [omesh] = build_omesh(geo)
% Construct an O-mesh consisting of 4 geometrical arc-segments surrounding a square
% center segment.
% The resulting Multi-Mesh object has the same properties as a single mesh
% object with the addition of the field "submeshes", which contain the
% individual meshes (in this case 5)

n_sub = 5; % Number of Patches to actually connect into one
% May be changed to only use a subset of the patches. 
% Attention: Expect translation in the solution in these cases.


% build iga mesh structure
num_meshes = min(length(geo), n_sub);
mesh = cell(1,num_meshes);
for i = 1:num_meshes
    mesh{1,i} = build_iga_mesh( geo{1,i} );   
    %mesh{1, i}.dim = 3;
    mesh{1,i}.gloElNodeCnt = mesh{1,i}.elNodeCnt;
    mesh{1,i}.nodeNum = zeros(mesh{1,i}.nCptsV, mesh{1,i}.nCptsU);
end

% assuming all circle segments have the same rows and columns
n_rows = mesh{1,1}.nCptsV; % Num Points along V-Axis
n_cols = mesh{1,1}.nCptsU; % Num Points along U-Axis
p = mesh{1,1}.p;
q = mesh{1,1}.q;
pq = (p+1)*(q+1);

% Overlap of ElCpts for number of connected elements
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
            %mesh{1,k}.nodeNum = zeros(n_rows, n_cols);
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
    mesh{1,5}.nodeNum(1, :) = mesh{1,1}.nodeNum(end, :); % First row of Pathc5 = last row of Patch1
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

%% GLobal Connectivity

global_node_id = 0;
global_ElNodeCnt = zeros(n_elems_total, (p+1)*(q+1));
global_elDoma = zeros(n_elems_total, 4);

for kk = 1:num_meshes
    % Speicher vorallokieren für die Konnektivitätsmatrix des aktuellen Patches
    mesh{1,kk}.gloElNodeCnt = zeros(mesh{1,kk}.nElems, mesh{1,kk}.nElemCpts);
    
    % Zähler für die Elemente des aktuellen Patches
    el_idx = 0;
    
    % Schleife über alle Elemente des aktuellen Patches
    for j = 1:mesh{1,kk}.nElemV
        for i = 1:mesh{1,kk}.nElemU
            global_node_id = global_node_id + 1; % Node index GLOBALLY
            el_idx = el_idx + 1; % Node index in THIS Mesh
            
            % Das Fenster greift auf die fertig gekoppelte nodeNum des Patches zu
            % Dank deiner Vorarbeit enthält diese bereits die globalen IDs!
            local_window = mesh{1,kk}.nodeNum(j:j+q, i:i+p);
            
            % Zeilenweise flachklopfen und in die globale Matrix schreiben
            mesh{1,kk}.gloElNodeCnt(el_idx, :) = reshape(local_window', 1, []);
            
            % TODO: Check if this works
            %mesh{1, kk}.elNodeCnt(el_idx, :) = reshape(local_window', 1, []);
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
%omesh.k = mesh{1,1}.p; % Copy to use for 3D Mesh Visualiation options
omesh.nCpts = n_points_total; % All available unique control points?
omesh.coords = globCoords;
omesh.initcoords = globCoords;
omesh.nElems = n_elems_total;
omesh.elDoma = global_elDoma;
omesh.nElemCpts = pq;
omesh.elNodeCnt = global_ElNodeCnt;
omesh.uKnots = mesh{1,1}.uKnots;
omesh.vKnots = mesh{1,1}.vKnots;
%omesh.wKnots = zeros(size(mesh{1,1}.vKnots)); % Copy to use for 3D Mesh Visualiation options
omesh.nCptsU = n_cols;
omesh.nCptsV = n_rows;
omesh.submeshes = mesh;

end