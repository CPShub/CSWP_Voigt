
% define the geometry
plate_circle = geo_circle_with_square([0,0], 1, 0.4);
%plate_circle = geo_circle_with_square_3d([0,0], 1, 0.4,0);

% Create the Mesh
multi_mesh = build_omesh(plate_circle);

dof = 3;
mat = default_mat();
mat.index = 114;
ndofs = dof * multi_mesh.nCpts;   % total dofs
u = zeros(ndofs + 6,1);
eps0 = [0,0.005,0]';
k0 = [0,0,0]';
curtime=1;

filename_pk2 = 'trash2';
fout = fopen(get_output_file_name(filename_pk2), 'w');
geo_s = geo_square();
mesh_s = build_iga_mesh(geo_s);

%[ Kglob, Rglob ] = globalstiffness_CSWP_PK2_meshwise(30, plate_circle, multi_mesh, mat, u , curtime,eps0,k0);
%[nliga_return] = nliga_returns( 30, geo_s, mesh_s, mat, [], [], fout,eps0, k0 );
[nliga_return] = nliga_returns( 30, plate_circle, multi_mesh, mat, [], [], fout,eps0, k0 );

%% interpolate the results on a grid along the X-Connection boundary

% Read the trash2 file
%vmesh = read_visual_mesh("trash2.msh");
u = nliga_return.u;
dof = 3;

% Compute DEF and STRESS results
boundary_elements = [13, 14, 15, 16, 65, 66, 67, 68];
M1_element_ids = [13, 14, 15, 16];
M5_element_ids = [65, 66, 67, 68];
boundary_elNodeCnt = multi_mesh.elNodeCnt(boundary_elements, :);
boundary_coords = multi_mesh.coords(boundary_elNodeCnt, :);


% Define Evaluateion grid in param space

% Mesh 1 Along x-axis from 0 to 1, y-axis along 0.9 (between 0.75 and 1)
n_eval_x = 20; % Num Points along X-Axis
n_eval_y = 10; % Num Points along Y-Axis
n_eval = n_eval_x * n_eval_y;

M1_range = [0.8, 1];
M5_range = [0, 0.2];


M1_evalpoints_param_x = linspace(0, 1, n_eval_x);
M1_evalpoints_param_y = linspace(M1_range(1), M1_range(2), n_eval_y);

doms = multi_mesh.elDoma(M1_element_ids, :);
is_inside = (M1_evalpoints_param_x >= doms(:,1)) & (M1_evalpoints_param_x < doms(:,2));
is_inside(end, M1_evalpoints_param_x == doms(end,2)) = true;
[row_idx, ~] = find(is_inside);
M1_evalpoints_e = M1_element_ids(row_idx);
M1_evalpoints_x = zeros(n_eval, 3);
M1_evalpoints_vms = zeros(n_eval, 1);
M1_evalpoints_ccy = zeros(n_eval, 6);
M1_evalpoints_dxalpha = zeros(n_eval, 6);
M1_evalpoints_ders = zeros(n_eval, 16);
M1_evalpoints_edsp = zeros(n_eval, 16);

% ==========


M5_evalpoints_param_x = linspace(0, 1, n_eval_x);
M5_evalpoints_param_y = linspace(M5_range(1), M5_range(2), n_eval_y);
doms = multi_mesh.elDoma(M5_element_ids, :);
is_inside = (M5_evalpoints_param_x >= doms(:,1)) & (M5_evalpoints_param_x < doms(:,2));
is_inside(end, M5_evalpoints_param_x == doms(end,2)) = true; % Randfall x=1.0 abfangen

% 3. Mappe die gefundenen Positionen zurück auf deine element_ids
[row_idx, ~] = find(is_inside);
M5_evalpoints_e = M5_element_ids(row_idx);
M5_evalpoints_x = zeros(n_eval, 3);
M5_evalpoints_vms = zeros(n_eval, 1);
M5_evalpoints_ccy = zeros(n_eval, 6);
M5_evalpoints_dxalpha = zeros(n_eval, 6);
M5_evalpoints_ders = zeros(n_eval, 16);
M5_evalpoints_edsp = zeros(n_eval, 16);


% Compute interpolation points for Mesh M1
sub_mesh = multi_mesh.submeshes{1,1};
j = 0;
for iy = 1:n_eval_y
    for ix = 1:n_eval_x
        p = [M1_evalpoints_param_x(ix), M1_evalpoints_param_y(iy)];
        e = M1_evalpoints_e(ix);
    
        sctr = multi_mesh.elNodeCnt(e,:);     % element control points index
        %exyz = sub_mesh.coords(sctr,:);  % element control points' coordinates
        nn = length(sctr);   % number of control points in the element
        nn3 = nn*3;          % degree of freedom of control points
        %nn3 = nn*2;
        elDoma = multi_mesh.elDoma(e, :);
        
        
        % Check if global globElNodeCnt should be used
        sctr = multi_mesh.elNodeCnt(e, :);
        sctrB = zeros(1, nn3);      
        sctrB(1:3:nn3) = 3*sctr - 2;% displacement in x direction
        sctrB(2:3:nn3) = 3*sctr-1;  % displacement in y direction
        sctrB(3:3:end) = 3*sctr;    % displacement in z direction
        
        edsp = u(sctrB);
        edsp = reshape(edsp, 3, nn);
    
        elCpts0 = sub_mesh.initcoords(sctr,:); % initial coordinates of el cont points
        elCpts(:,1:3)=elCpts0(:,1:3)+edsp';
    
    
        % Evalute in the interpolation grid
        [N,ders] = nurbs_derivatives(p, plate_circle, sub_mesh);
        jmatrix = ders*elCpts0(:,1:dof-1); %Because the mapping is in 2D
        ders =  jmatrix \ ders;    
        ders3D = zeros(3,size(elCpts,1));
        ders3D(1:2,:) = ders;
        x = N.*elCpts(:,1:dof)';
        x = sum(x,2);
    
        dx_alpha = edsp * ders3D';
        F = def_gradient(eps0, k0, x, dx_alpha);
    
        % Retrieve PK2 material response as PK2 Stress and dtangent 
        [ stress, ~ ] = material_CSWP_PK2_hyperelasticity( dof, mat, F );
        
        ccy = pk2cauchy(stress, F);
            
        % Store the interpolated positional and stress values
        j = j + 1;
        M1_evalpoints_x(j, :) = x;
        M1_evalpoints_vms(j, 1) = von_mises(ccy');
        M1_evalpoints_ccy(j, :) = ccy;
        M1_evalpoints_dxalpha(j, :) = reshape(dx_alpha(:, 1:2), 1, 6);
        M1_evalpoints_ders(j, :) = ders(2, :);
        M1_evalpoints_edsp(j, :) = edsp(3, :);
    end
end

% Do the same for M5
sub_mesh = multi_mesh.submeshes{1,5};
j = 0;
for iy = 1:n_eval_y
    for ix = 1:n_eval_x
        p = [M5_evalpoints_param_x(ix), M5_evalpoints_param_y(iy)];
        e = M5_evalpoints_e(ix);
    
        sctr = multi_mesh.elNodeCnt(e,:);     % element control points index
        %exyz = sub_mesh.coords(sctr,:);  % element control points' coordinates
        nn = length(sctr);   % number of control points in the element
        nn3 = nn*3;          % degree of freedom of control points
        %nn3 = nn*2;
        elDoma = multi_mesh.elDoma(e, :);
        
        
        % Check if global globElNodeCnt should be used
        sctr = multi_mesh.elNodeCnt(e, :);
        sctrB = zeros(1, nn3);      
        sctrB(1:3:nn3) = 3*sctr - 2;% displacement in x direction
        sctrB(2:3:nn3) = 3*sctr-1;  % displacement in y direction
        sctrB(3:3:end) = 3*sctr;    % displacement in z direction
        
        edsp = u(sctrB);
        edsp = reshape(edsp, 3, nn);
    
        elCpts0 = multi_mesh.initcoords(sctr,:); % initial coordinates of el cont points
        elCpts(:,1:3)=elCpts0(:,1:3)+edsp';
    
    
        % Evalute in the interpolation grid
        [N,ders] = nurbs_derivatives(p, plate_circle, sub_mesh);
        jmatrix = ders*elCpts0(:,1:dof-1); %Because the mapping is in 2D
        ders =  jmatrix \ ders;    
        ders3D = zeros(3,size(elCpts,1));
        ders3D(1:2,:) = ders;
        x = N.*elCpts(:,1:dof)';
        x = sum(x,2);
    
        dx_alpha = edsp * ders3D';
        F = def_gradient(eps0, k0, x, dx_alpha);
    
        % Retrieve PK2 material response as PK2 Stress and dtangent 
        [ stress, ~ ] = material_CSWP_PK2_hyperelasticity( dof, mat, F );
        
        ccy = pk2cauchy(stress, F);
            
        % Store the interpolated positional and stress values
        j = j + 1;
        M5_evalpoints_x(j, :) = x;
        M5_evalpoints_vms(j, 1) = von_mises(ccy');
        M5_evalpoints_ccy(j, :) = ccy;
        M5_evalpoints_dxalpha(j, :) = reshape(dx_alpha(:, 1:2), 1, 6);
        M5_evalpoints_ders(j, :) = ders(2, :);
        M5_evalpoints_edsp(j, :) = edsp(3, :);
    end
end

% Visualize
disp("")

figure
grid on
hold on

% plot the Controlpoints
scatter(multi_mesh.coords(multi_mesh.elNodeCnt([13, 14, 15, 16], :), 1), multi_mesh.coords(multi_mesh.elNodeCnt([13, 14, 15, 16], :), 2), 20, "red", "filled");
scatter(multi_mesh.coords(multi_mesh.elNodeCnt([65, 66, 67, 68], :), 1), multi_mesh.coords(multi_mesh.elNodeCnt([65, 66, 67, 68], :), 2), 20, "blue", "filled");

% plot the interpolation points in between
scatter(M1_evalpoints_x(:, 1), M1_evalpoints_x(:, 2), 20, "green", "filled");



% Figure in 3D
% figure
% 
% % plot the Controlpoints
% scatter3(multi_mesh.coords(multi_mesh.elNodeCnt([13, 14, 15, 16], :), 1), ...
%     multi_mesh.coords(multi_mesh.elNodeCnt([13, 14, 15, 16], :), 2), ...
%     zeros(25, 1), ...
%     20, "red", "filled");
% grid on
% hold on
% scatter3(multi_mesh.coords(multi_mesh.elNodeCnt([65, 66, 67, 68], :), 1), ...
%     multi_mesh.coords(multi_mesh.elNodeCnt([65, 66, 67, 68], :), 2), ...
%     zeros(25, 1), ...
%     20, "blue", "filled");



M1_last_ids = 181:200; % TODO: Make Automatic
M5_first_ids = 1:20; % TODO: Make Automatic


figure
subtitle("Mesh 1 dx_alpha6 Component analysis")
subplot(2,2,2)
hold on
grid on
imagesc(M1_evalpoints_ders(M1_last_ids, 9:end)')
colorbar()
yyaxis right
plot((1:20)', sum(M1_evalpoints_ders(M1_last_ids, 9:end)'))
subplot(2,2,3)
hold on
grid on
imagesc(M1_evalpoints_edsp(M1_last_ids, 9:end))
colorbar()
subplot(2,2,4)
hold on
grid on
M1_dxalpha_matrix = flip(diag(M1_evalpoints_dxalpha(M1_last_ids, 6)));
imagesc(M1_dxalpha_matrix);
colorbar()
clim([min(M1_evalpoints_dxalpha(M1_last_ids, 6)), max(M1_evalpoints_dxalpha(M1_last_ids, 6))])
%scatter((1:20'), M1_evalpoints_dxalpha(M1_last_ids, 6), 20, 'MarkerEdgeColor', "red", "DisplayName", sprintf("Dx_alpha(6)"));


figure
subtitle("Mesh 5 dx_alpha6 Component analysis")
subplot(2,2,2)
hold on
grid on
imagesc(M5_evalpoints_ders(M5_first_ids, :)')
colorbar()
yyaxis right
plot(M5_first_ids', sum(M1_evalpoints_ders(M5_first_ids, 9:end)'))
subplot(2,2,3)
hold on
grid on
imagesc(M5_evalpoints_edsp(M5_first_ids, :))
colorbar()
subplot(2,2,4)
hold on
grid on
M5_dxalpha_matrix = flip(diag(M5_evalpoints_dxalpha(M5_first_ids, 6)));
imagesc(M5_dxalpha_matrix);
colorbar()
clim([min(M5_evalpoints_dxalpha(M5_first_ids, 6)), max(M5_evalpoints_dxalpha(M5_first_ids, 6))])
%scatter((1:20'), M1_evalpoints_dxalpha(M1_last_ids, 6), 20, 'MarkerEdgeColor', "red", "DisplayName", sprintf("Dx_alpha(6)"));











%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DX_ALPHA
figure
M1_last_ids = 181:200; % TODO: Make Automatic
M5_first_ids = 1:20; % TODO: Make Automatic

subplot(2,1,1);
grid on
hold on
colors = ["#C44E52", "#55A868", "#4C72B0", "#DD8452", "#8172B3", "#937860"];
rgb_matrix = validatecolor(colors, 'multiple');

title("M1 - Dx_alpha")
for kk = 1:6
    % Iterate over all ccy components
    scatter(M1_evalpoints_param_x, M1_evalpoints_dxalpha(M1_last_ids, kk), ...
        20, 'MarkerEdgeColor', rgb_matrix(kk,:), ...
        "DisplayName", sprintf("Dx_alpha(%d)", kk));
    plot(M1_evalpoints_param_x, M1_evalpoints_dxalpha(M1_last_ids, kk), "Color", rgb_matrix(kk,:))
end

subplot(2,1,2)
title("M5")
grid on
hold on
for kk = 1:6
    % Iterate over all ccy components
    scatter(M5_evalpoints_param_x, M5_evalpoints_dxalpha(M5_first_ids, kk), ...
        20, 'MarkerEdgeColor', rgb_matrix(kk,:), ...
        "DisplayName", sprintf("Dx_alpha: %d", kk));
    plot(M5_evalpoints_param_x, M5_evalpoints_dxalpha(M5_first_ids, kk), ...
        "Color", rgb_matrix(kk,:))
end
legend();


























%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CCY
figure
%suptitle("CCY Components over Mesh Boundary")
% only plot the very Last of M1 and very First of M5 in a 2D plot

M1_last_ids = 181:200; % TODO: Make Automatic
M5_first_ids = 1:20; % TODO: Make Automatic

subplot(2,1,1);
grid on
hold on
colors = ["#C44E52", "#55A868", "#4C72B0", "#DD8452", "#8172B3", "#937860"];
rgb_matrix = validatecolor(colors, 'multiple');

title("M1")
for kk = 5:5
    % Iterate over all ccy components
    scatter(M1_evalpoints_param_x, M1_evalpoints_ccy(M1_last_ids, kk), ...
        20, 'MarkerEdgeColor', rgb_matrix(kk,:), ...
        "DisplayName", sprintf("CCY: %d", kk));
    plot(M1_evalpoints_param_x, M1_evalpoints_ccy(M1_last_ids, kk), "Color", rgb_matrix(kk,:))
end

subplot(2,1,2)
title("M5")
grid on
hold on
for kk = 5:5
    % Iterate over all ccy components
    scatter(M5_evalpoints_param_x, M5_evalpoints_ccy(M5_first_ids, kk), ...
        20, 'MarkerEdgeColor', rgb_matrix(kk,:), ...
        "DisplayName", sprintf("CCY: %d", kk));
    plot(M5_evalpoints_param_x, M5_evalpoints_ccy(M5_first_ids, kk), ...
        "Color", rgb_matrix(kk,:))
end
legend();


%%%%%%%%%%%%%%%%%%
% POSITION


figure
title("Position over Mesh Boundary")
scatter3( ...
    M1_evalpoints_x(:, 1), ...
    M1_evalpoints_x(:, 2), ...
    M1_evalpoints_x(:, 3), ...
    20, "green", "filled");
grid on
hold on
scatter3( ...
    M5_evalpoints_x(:, 1), ...
    M5_evalpoints_x(:, 2), ...
    M5_evalpoints_x(:, 3), ...
    50, "magenta");
xlabel("X")
ylabel("Y")








% plot the interpolation points in between
figure
title("Von-Mises Stress over the Mesh Boundary")
scatter3( ...
    M1_evalpoints_x(:, 1), ...
    M1_evalpoints_x(:, 2), ...
    M1_evalpoints_vms, ...
    20, "green", "filled");
grid on
hold on
scatter3( ...
    M5_evalpoints_x(:, 1), ...
    M5_evalpoints_x(:, 2), ...
    M5_evalpoints_vms, ...
    20, "magenta", "filled");
xlabel("X")
ylabel("Y")


