
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
%eps0 = [0.005,0,0]'; 
eps0 = [0,0.01,0]';
%eps0 = [0,0,0]';
k0 = [0,0,0]';
%k0 = [0, 0, 0.05]';
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
M2_element_ids = [17, 21, 25, 29];


% Define Evaluateion grid in param space

% Mesh 1 Along x-axis from 0 to 1, y-axis along 0.9 (between 0.75 and 1)
n_eval_x = 20; % Num Points along X-Axis
n_eval_y = 10; % Num Points along Y-Axis
n_eval = n_eval_x * n_eval_y;

M1_range = [0.8, 1];
M5_range = [0, 0.2];
M2_range = [0, 0.2];


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

M1_evalpoints_dXdx1s = zeros(n_eval, 3);
M1_evalpoints_dXdx2s = zeros(n_eval, 3);



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


% ===========


M2_evalpoints_param_x = linspace(M2_range(1), M2_range(2), n_eval_y);
M2_evalpoints_param_y = linspace(0, 1, n_eval_x);
doms = multi_mesh.elDoma(M2_element_ids, :);
is_inside = (M2_evalpoints_param_y >= doms(:,3)) & (M2_evalpoints_param_y < doms(:,4));
is_inside(end, end) = true; % Randfall y=1.0 abfangen

% 3. Mappe die gefundenen Positionen zurück auf deine element_ids
[row_idx, ~] = find(is_inside);
M2_evalpoints_e = M2_element_ids(row_idx);
M2_evalpoints_x = zeros(n_eval, 3);
M2_evalpoints_vms = zeros(n_eval, 1);
M2_evalpoints_ccy = zeros(n_eval, 6);
M2_evalpoints_dxalpha = zeros(n_eval, 6);
M2_evalpoints_ders = zeros(n_eval, 16);
M2_evalpoints_edsp = zeros(n_eval, 16);




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
    
        dXdx1 = sum(ders(1, :) .* elCpts(:, 1:dof)', 2);
        dXdx2 = sum(ders(2, :) .* elCpts(:, 1:dof)', 2);

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

        M1_evalpoints_dXdx1s(j, :) = dXdx1;
        M1_evalpoints_dXdx2s(j, :) = dXdx2;
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

% Do the same again for Mesh 2
sub_mesh = multi_mesh.submeshes{1,2};
j = 0;
for iy = 1:n_eval_x
    for ix = 1:n_eval_y
        p = [M2_evalpoints_param_x(ix), M2_evalpoints_param_y(iy)];
        e = M2_evalpoints_e(iy);

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
        M2_evalpoints_x(j, :) = x;
        M2_evalpoints_vms(j, 1) = von_mises(ccy');
        M2_evalpoints_ccy(j, :) = ccy;
        M2_evalpoints_dxalpha(j, :) = reshape(dx_alpha(:, 1:2), 1, 6);
        M2_evalpoints_ders(j, :) = ders(2, :);
        M2_evalpoints_edsp(j, :) = edsp(3, :);
    end
end



% Visualize
% disp("")
% 
% figure
% grid on
% hold on
% axis square
% 
% % plot the Controlpoints
% scatter(multi_mesh.coords(multi_mesh.elNodeCnt(M1_element_ids, :), 1), multi_mesh.coords(multi_mesh.elNodeCnt(M1_element_ids, :), 2), 20, "red", "filled");
% scatter(multi_mesh.coords(multi_mesh.elNodeCnt(M5_element_ids, :), 1), multi_mesh.coords(multi_mesh.elNodeCnt(M5_element_ids, :), 2), 20, "blue", "filled");
% scatter(multi_mesh.coords(multi_mesh.elNodeCnt(M2_element_ids, :), 1), multi_mesh.coords(multi_mesh.elNodeCnt(M2_element_ids, :), 2), 20, "magenta", "filled");
% 
% % plot the interpolation points in between
% scatter(M1_evalpoints_x(:, 1), M1_evalpoints_x(:, 2), 20, "green", "filled");
% scatter(M5_evalpoints_x(:, 1), M5_evalpoints_x(:, 2), 20, "black", "filled");
% scatter(M2_evalpoints_x(:, 1), M2_evalpoints_x(:, 2), 20, "cyan", "filled");
% 
% 
% 
% % Quiver Plot the of the local derivatives of Position
% % Currently only for element 1
% figure
% 
% scatter3(M1_evalpoints_x(:, 1), M1_evalpoints_x(:, 2), M1_evalpoints_x(:, 3), 20, "green", "filled")
% grid on
% hold on
% quiver3(M1_evalpoints_x(:, 1), M1_evalpoints_x(:, 2), M1_evalpoints_x(:, 3, ...
% M1_evalpoints_dXdx1s(:, 1), M1_evalpoints_dXdx1s(:, 2), M1_evalpoints_dXdx1s(:, 3), ...
% "red"))
% quiver3(M1_evalpoints_x(:, 1), M1_evalpoints_x(:, 2), M1_evalpoints_x(:, 3, ...
%     M1_evalpoints_dXdx2s(:, 1), M1_evalpoints_dXdx2s(:, 2), M1_evalpoints_dXdx2s(:, 3), ...
%     "blue"))
% 
% xlabel("X")
% ylabel("Y")









%% 3D Correlation Analysis between edsp and ders



% Plot the Controlpoint positions for the first n=1:20 visualization-nodes
% in M1
M1_last_ids = 181:200; 
M5_first_ids = 1:20; 

fig = figure;
fig.Position = [100, 100, 800, 600];

ax = axes('Parent', fig, 'Position', [0.1, 0.25, 0.8, 0.68]);

sld = uicontrol('Parent', fig, 'Style', 'slider', ...
    'Min', 1, 'Max', 20, 'Value', 1, ...
    'SliderStep', [1/19, 1/19], ...
    'Position', [150, 40, 500, 20]);

txt = uicontrol('Parent', fig, 'Style', 'text', ...
    'Position', [350, 70, 100, 20], ...
    'String', 'Index i = 1');

update_plot(1, ax, multi_mesh, M1_evalpoints_e, M1_evalpoints_edsp, M1_evalpoints_ders, M1_evalpoints_x, M1_evalpoints_dxalpha);

sld.Callback = @(src, event) callback_function(src, txt, ax, multi_mesh, M1_evalpoints_e, M1_evalpoints_edsp, M1_evalpoints_ders, M1_evalpoints_x, M1_evalpoints_dxalpha);

function callback_function(src, txt_handle, ax_handle, multi_mesh, M1_evalpoints_e, M1_evalpoints_edsp, M1_evalpoints_ders, M1_evalpoints_x, M1_evalpoints_dxalpha)
    i = round(src.Value);
    src.Value = i;
    txt_handle.String = ['Index i = ', num2str(i)];
    update_plot(i, ax_handle, multi_mesh, M1_evalpoints_e, M1_evalpoints_edsp, M1_evalpoints_ders, M1_evalpoints_x, M1_evalpoints_dxalpha);
end

function update_plot(i, ax_handle, multi_mesh, M1_evalpoints_e, M1_evalpoints_edsp, M1_evalpoints_ders, M1_evalpoints_x, M1_evalpoints_dxalpha)
    cla(ax_handle);
    hold(ax_handle, 'on');
    grid(ax_handle, 'on');
    
    el_i = M1_evalpoints_e(i);
    i = i + 180;
    cnt_x = multi_mesh.initcoords(multi_mesh.elNodeCnt(el_i, :), :);
    edsp3_i = M1_evalpoints_edsp(i, :);
    ders2_i = M1_evalpoints_ders(i, :);
    dx_alpha6_i = edsp3_i .* ders2_i;
    x_i = M1_evalpoints_x(i, :);
    dx_alpha6_INT_i = M1_evalpoints_dxalpha(i, 6);
    
    scale1 = max(abs(ders2_i)) / max(abs(edsp3_i));
    scale2 = max(abs(ders2_i)) / max(abs(dx_alpha6_i));
    scale3 = max(abs(ders2_i)) / abs(x_i(3));
    scale4 = max(abs(ders2_i)) / abs(dx_alpha6_INT_i);
    
    z_lower_cnt = scale1 * edsp3_i;
    z_upper_cnt = scale2 * dx_alpha6_i;
    
    for jj = 1:16
        plot3(ax_handle, [cnt_x(jj, 1), cnt_x(jj, 1)], [cnt_x(jj, 2), cnt_x(jj, 2)], [z_lower_cnt(jj), z_upper_cnt(jj)], "--", "Color", [0.3, 0.3, 0.3], 'HandleVisibility','off');
    end
    
    z_lower_int = scale3 * x_i(:, 3);
    z_upper_int = scale4 * dx_alpha6_INT_i;
    plot3(ax_handle, [x_i(:, 1), x_i(:, 1)], [x_i(:, 2), x_i(:, 2)], [z_lower_int, z_upper_int], "-", "Color", "black", "LineWidth", 2, 'HandleVisibility','off');
    
    scatter3(ax_handle, cnt_x(:, 1), cnt_x(:, 2), scale1 * edsp3_i, 30, "red", "filled", "^", "DisplayName", "CNT Z-Displacement"); 
    scatter3(ax_handle, cnt_x(:, 1), cnt_x(:, 2), ders2_i, 60, "blue", "o", "DisplayName", "CNT Y-Ders"); 
    scatter3(ax_handle, cnt_x(:, 1), cnt_x(:, 2), scale2 * dx_alpha6_i, 30, "green", "filled", "square", "DisplayName", "CNT DX_alpha(6)"); 
    scatter3(ax_handle, x_i(:, 1), x_i(:, 2), scale3 * x_i(:, 3), 30, "black", "filled", "*", 'DisplayName', "INT Position");
    scatter3(ax_handle, x_i(:, 1), x_i(:, 2), scale4 * dx_alpha6_INT_i, 30, "magenta", "square", 'DisplayName', "INT DX_alpha(6)");
    
    title(ax_handle, ["Various Variables over Associated Control Points for i = ", num2str(i)]);
    legend(ax_handle, 'show');
    xlabel(ax_handle, "X");
    ylabel(ax_handle, "Y");
    hold(ax_handle, 'off');
end