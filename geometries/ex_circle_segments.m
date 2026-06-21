% define the geometry
plate_circle = geo_circle_with_square([0,0], 1, 0.4);
n_sub = 5; % Number of Patches to actually connect into one


% build iga mesh structure
num_meshes = min(length(plate_circle), n_sub);
mesh = cell(1,num_meshes);
for i = 1:num_meshes
    mesh{1,i} = build_iga_mesh( plate_circle{1,i} );   
    mesh{1,i}.gloElNodeCnt = mesh{1,i}.elNodeCnt;
    mesh{1,i}.nodeNum = zeros(mesh{1,i}.nCptsV, mesh{1,i}.nCptsU);
end

% assuming all circle segments have the same rows and columns
n_rows = mesh{1,1}.nCptsV;
n_cols = mesh{1,1}.nCptsU;
p = mesh{1,1}.p;
q = mesh{1,1}.q;
pq = (p+1)*(q+1);

% Overlap of ElCpts for number of connected elements
rem = [0, 
    mesh{1,1}.nCptsV,
    2*mesh{1,1}.nCptsV,
    4*mesh{1,1}.nCptsV,
    4*mesh{1,1}.nCptsV + 4*(mesh{1,1}.nCptsU-1)];

n_points_total = num_meshes * (mesh{1,1}.nCpts) - rem(num_meshes);
%n_points_total = (n_rows * (n_cols-1))*4 + (n_cols-2)*(n_rows-2); % Total of unique points
n_elems_total = num_meshes * mesh{1,1}.nElems;

% node numbering (top left start, row-wise)
for k = 1:num_meshes
    for j = 1:n_rows % Radius Count
        for i = 1:n_cols % Radial Count
            mesh{1,k}.nodeNum(j,i) = (j-1)*(mesh{1,k}.nCptsU)+i;
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
    for j = 2:mesh{1,2}.nCptsV % Rows going up (+Y)
        for i = 1:mesh{1,2}.nCptsU % Columns going right (+X)
            bb = (j-1)*(mesh{1,2}.nCptsV)+i;    % Local Node Index
            nodal_num = nodal_num + 1;          % global Node Index
            
            mesh{1,2}.nodeNum(j,i) = nodal_num; % Assign correct GLOBAL Node Index
            globCoords(nodal_num, :) = mesh{1,2}.coords(bb, :); % Assign Local coords as global coords (Index-Shift)
        end
    end
end

if num_meshes >= 3
    % node Numbering & globCoords assignment for patch 3
    mesh{1,3}.nodeNum(:, end) = mesh{1,2}.nodeNum(end, :); % Last Row of Patch2 = Last Column of Patch3
    for j = 1:mesh{1,3}.nCptsV % Rows going up (+Y)
        for i = 1:mesh{1,3}.nCptsU-1 % Columns going right (+X)
            bb = (j-1)*(mesh{1,3}.nCptsV)+i;    % Local Node Index
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
    for j = 2:mesh{1,4}.nCptsV-1 % Rows going up (+Y)
        for i = 1:mesh{1,4}.nCptsU % Columns going right (+X)
            bb = (j-1)*(mesh{1,4}.nCptsV)+i;    % Local Node Index
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
    
    for j = 2:mesh{1,4}.nCptsV-1 % Rows going up (+Y)
        for i = 2:mesh{1,4}.nCptsU-1 % Columns going right (+X)
            bb = (j-1)*(mesh{1,5}.nCptsV)+i;    % Local Node Index
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
            global_ElNodeCnt(global_node_id, :) = reshape(local_window', 1, []);
            global_elDoma(global_node_id, :) = mesh{1, kk}.elDoma(el_idx, :);
            
        end
    end
end


%% Finalize into Mesh form
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


multi_mesh.dim = mesh{1,1}.dim;
multi_mesh.p = mesh{1,1}.p;
multi_mesh.q = mesh{1,1}.q;

multi_mesh.nCpts = n_points_total; % All available unique control points?
multi_mesh.coords = globCoords;
multi_mesh.initcoords = globCoords;
multi_mesh.nElems = n_elems_total;
multi_mesh.elDoma = global_elDoma;
multi_mesh.nElemCpts = pq;
multi_mesh.elNodeCnt = global_ElNodeCnt;
multi_mesh.uKnots = mesh{1,1}.uKnots;
multi_mesh.vKnots = mesh{1,1}.vKnots;
multi_mesh.nCptsU = mesh{1,1}.nCptsU;
multi_mesh.nCptsV = mesh{1,1}.nCptsV;
multi_mesh.submeshes = mesh;

multi_geo = cell(1, num_meshes);
%multi_geo.form = "B-NURBS";
%multi_geo.dim = 4;
for jjj = 1:num_meshes
    multi_geo{1, jjj} = plate_circle{1, jjj};
end
%multi_mesh.submesh_node_ids = global_ElNodeCnt.reshape(mesh{1,1}.nElems, 5, :);% All global node Indicees associated with Meshes 1 to 5 as separate lists
%mutli_mesh.submesh_element_ids = ;% All global Element Indicees associated ...





%% Compute a globalstiffness-matrix using this new multi-mesh

new_school = 1;

if new_school
    dof = 3;
    mat = default_mat();
    mat.index = 114;
    ndofs = dof * multi_mesh.nCpts;   % total dofs
    u = zeros(ndofs + 6,1);
    eps0 = [0,0,0]';
    k0 = [0,0,0.05]';
    curtime=1;

    filename_pk2 = 'trash';
    fout = fopen(get_output_file_name(filename_pk2), 'w');
    geo_s = geo_square();
    mesh_s = build_iga_mesh(geo_s);

    %[ Kglob, Rglob ] = globalstiffness_CSWP_PK2_meshwise(30, plate_circle, multi_mesh, mat, u , curtime,eps0,k0);
    %[nliga_return] = nliga_returns( 30, geo_s, mesh_s, mat, [], [], fout,eps0, k0 );
    [nliga_return] = nliga_returns( 30, multi_geo, multi_mesh, mat, [], [], fout,eps0, k0 );
    %[nliga_return] = nliga_returns(30, plate_circle{1,1}, mesh{1,1}, mat, [], [], fout, eps0, k0);
    %[nliga_return] = nliga_returns(30, geo_sub, mesh_sub, mat, [], [], fout, eps0, k0); 
    
    
    %% Visualization
    plot_u_over_time(multi_mesh, {}, nliga_return.u, "", {});

    % Nodes
    figure('Color', 'w');
    hold on;
    axis equal;
    grid on;
    set(gca, 'FontName', 'Helvetica', 'FontSize', 11);
    scatter(multi_mesh.coords(:, 1), multi_mesh.coords(:, 2), 30, "blue", 'filled');
    for rr = 1:multi_mesh.nCpts
        text(multi_mesh.coords(rr, 1), multi_mesh.coords(rr, 2), sprintf(' %d', rr), ...
            'Color', "blue", 'FontSize', 9, 'VerticalAlignment', 'bottom');
    end
    
end