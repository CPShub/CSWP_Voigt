% define the geometry
plate_circle = geo_circle_with_square([0,0], 1, 0.4);


% build iga mesh structure
mesh = cell(1,length(plate_circle));
for i = 1:length(plate_circle)
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

n_points_total = (n_rows * (n_cols-1))*4 + (n_cols-2)*(n_rows-2); % Total of unique points
n_elems_total = 5 * mesh{1,1}.nElems;

% node numbering (top left start, row-wise)
for k = 1:5
    for j = 1:n_rows % Radius Count
        for i = 1:n_cols % Radial Count
            mesh{1,k}.nodeNum(j,i) = (j-1)*(mesh{1,k}.nCptsU)+i;
        end
    end
end

% assign global coordinates of the first patch
%mesh{1,1}.globCoords = mesh{1,1}.coords; % Create 49 Global coords
globCoords = zeros(n_points_total, 4);
globCoords(1:(n_rows*n_cols) , :) = mesh{1,1}.coords;

%  node numbering for patch 2 to 3
max_num = n_rows*n_cols;
for kk = 2:3 
    mesh{1,kk}.nodeNum(:,1) = mesh{1,kk-1}.nodeNum(:,end);
    mesh{1,kk}.globCoords = zeros(mesh{1,kk}.nCptsV*(mesh{1,kk}.nCptsU-1),4); % Create an ADDITIONAL 42 unique Global Coords
    for j = 1:mesh{1,kk}.nCptsV
        for i = 2:mesh{1,kk}.nCptsU
            bb = (j-1)*(mesh{1,kk}.nCptsU-1)+i-1; % Local Node Index
            nodal_num = bb + max_num;%mesh{1,kk-1}.nodeNum(end,end); % global Node Index
            mesh{1,kk}.nodeNum(j,i) = nodal_num; % Assign correct GLOBAL Node Index
            globCoords(nodal_num, :) = mesh{1,kk}.coords(bb, :); % Assign Local coords as global coords (Index-Shift)
        end
    end
end


% Create only 35 unique global coords for patch 4
kk = 4;
mesh{1,kk}.nodeNum(:,1) = mesh{1,kk-1}.nodeNum(:,end);
bb = 0;
for j = 1:mesh{1,kk}.nCptsV
    for i = 2:mesh{1,kk}.nCptsU-1 % Not the First or last column
        bb = bb + 1; % Local Node Index
        nodal_num = bb + mesh{1,kk-1}.nodeNum(end,end); % global Node Index
        mesh{1,kk}.nodeNum(j,i) = nodal_num; % Assign correct GLOBAL Node Index
        globCoords(nodal_num, :) = mesh{1,kk}.coords(bb, :); % Assign Local coords as global coords (Index-Shift)
    end
end

%overwrite the last column with the node indecees from the first pathc
% This does NOT affect the global coords, because these are only registered
% for UNIQUE Points
mesh{1,4}.nodeNum(:, end) = mesh{1,1}.nodeNum(:, 1);


% Enumerate the nodes in Patch 5
bb = 0;
kk = 5;
for j = 2:mesh{1,5}.nCptsV-1 % For all inner rows
    for i = 2:mesh{1,5}.nCptsU-1 % For all inner columns
        bb = (j-1)*(mesh{1,kk}.nCptsU)+i;               % Local Node Index
        nodal_num = nodal_num + 1;   % Global Node Index

        mesh{1,kk}.nodeNum(j,i) = nodal_num; % Assign correct GLOBAL Node Index
        globCoords(nodal_num, :) = mesh{1,kk}.coords(bb, :); % Assign Local coords as global coords (Index-Shift)
    end
end

% Connect bottom(first) row of Patch 5 with top row of Patch 1
mesh{1,5}.nodeNum(1,:) = mesh{1,1}.nodeNum(1,:);

% Connect right(last) column of Patch 5 with left (first) row of Patch 1
mesh{1,5}.nodeNum(:,end) = mesh{1,2}.nodeNum(1,:)';

% Connect top(last) row of Patch 5 with first row of Patch 1
mesh{1,5}.nodeNum(end,:) = flip(mesh{1,3}.nodeNum(1, :));

mesh{1,5}.nodeNum(:,1) = flip(mesh{1,4}.nodeNum(1, :)');

%% GLobal Connectivity

global_node_id = 0;
global_ElNodeCnt = zeros(n_elems_total, (p+1)*(q+1));
global_elDoma = zeros(n_elems_total, 4);
for kk = 1:length(mesh)
    % Speicher vorallokieren für die Konnektivitätsmatrix des aktuellen Patches
    mesh{1,kk}.gloElNodeCnt = zeros(mesh{1,kk}.nElems, mesh{1,kk}.nElemCpts);
    
    % Zähler für die Elemente des aktuellen Patches
    el_idx = 0;
    
    % Schleife über alle Elemente des aktuellen Patches
    for j = 1:mesh{1,kk}.nElemV
        for i = 1:mesh{1,kk}.nElemU
            global_node_id = global_node_id + 1;
            el_idx = el_idx + 1;
            
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
multi_mesh.nElems = 5 * mesh{1,1}.nElems;
multi_mesh.elDoma = global_elDoma;
multi_mesh.nElemCpts = pq;
multi_mesh.elNodeCnt = global_ElNodeCnt;
multi_mesh.uKnots = mesh{1,1}.uKnots;
multi_mesh.vKnots = mesh{1,1}.vKnots;
multi_mesh.nCptsU = mesh{1,1}.nCptsU;
multi_mesh.nCptsV = mesh{1,1}.nCptsV;
%multi_mesh.submesh_node_ids = global_ElNodeCnt.reshape(mesh{1,1}.nElems, 5, :);% All global node Indicees associated with Meshes 1 to 5 as separate lists
%mutli_mesh.submesh_element_ids = ;% All global Element Indicees associated ...


%% Visualize as a checkup
figure('Color', 'w');
hold on;
axis equal;
grid on;
set(gca, 'FontName', 'Helvetica', 'FontSize', 11);
colors = ["blue", "red", "green", "violet", "black"];
offsets = [[0, -1]; [1,0];[0,1];[-1,0];[0,0]] .* 0.1;
for submesh_id = 1:5
    submesh_ids = (submesh_id-1)*49+1:(submesh_id)*49;
    submesh_coords = multi_mesh.coords(submesh_ids,1:2) + offsets(submesh_id, :);
    scatter(submesh_coords(:, 1), submesh_coords(:, 2), 30, colors(submesh_id), 'filled');
    for rr = 1:49
        n = submesh_ids(rr);
        text(submesh_coords(rr, 1), submesh_coords(rr, 2), sprintf(' %d', n), ...
            'Color', colors(submesh_id), 'FontSize', 9, 'VerticalAlignment', 'bottom');
    end
end
xlabel('X', 'FontWeight', 'bold');
ylabel('Y', 'FontWeight', 'bold');
title('Global Multi-Patch Mesh Verification', 'FontSize', 13, 'FontWeight', 'bold');
hold off;






figure('Color', 'w', 'Position', [100, 100, 900, 800]);
hold on;
axis equal;
grid on;
set(gca, 'FontName', 'Helvetica', 'FontSize', 11);

for e = 1:multi_mesh.nElems
    nodes = multi_mesh.elNodeCnt(e, :);
    elem_coords = multi_mesh.coords(nodes, 1:2);
    
    k = convhull(elem_coords(:,1), elem_coords(:,2));
    
    fill(elem_coords(k,1), elem_coords(k,2), rand(1,3), ...
        'FaceAlpha', 0.3, 'EdgeColor', [0.3 0.3 0.3], 'LineWidth', 1.0);
    
    cx = mean(elem_coords(:,1));
    cy = mean(elem_coords(:,2));
    
    text(cx, cy, num2str(e), 'Color', 'k', 'FontWeight', 'bold', ...
        'FontSize', 10, 'HorizontalAlignment', 'center');
end

xlabel('X', 'FontWeight', 'bold');
ylabel('Y', 'FontWeight', 'bold');
title('Global Multi-Patch Mesh Verification', 'FontSize', 13, 'FontWeight', 'bold');
hold off;














% Scatter the Nodes with associated global node Numbers as text





% draw the Elements as rectangles over their domains in random colors
% Each element has its element ID as a black bold text in the center


figure
hold on
grid on
i = 3;
labels = reshape(string(mesh{1, i}.nodeNum)', 49, 1);
scatter(mesh{1,i}.coords(:, 1), mesh{1,i}.coords(:, 2), 20, "red", "filled")
textscatter(mesh{1,i}.coords(:, 1), mesh{1,i}.coords(:, 2), labels)


figure
axis equal
axis square
hold on
grid on
colors = ["red", "green", "blue", "magenta", "black"];
for i = 1:5
    labels = reshape(string(mesh{1, i}.nodeNum)', 49, 1);
    scatter(mesh{1,i}.globCoords(:, 1), mesh{1,i}.globCoords(:, 2), 20, colors(i), "filled")
    textscatter(mesh{1,i}.globCoords(:, 1), mesh{1,i}.globCoords(:, 2), labels)
end