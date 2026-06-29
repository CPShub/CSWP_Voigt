
% Create the Mesh
geo = geo_circle_with_square([0,0], 1, 0.4);
multi_mesh = build_omesh(geo);

% Solve the Simulation and save to a file
dof = 3;
mat = default_mat();
mat.index = 114;
ndofs = dof * multi_mesh.nCpts;   % total dofs
u = zeros(ndofs + 6,1);
eps0 = [0,0.005,0]';
k0 = [0,0,0]';
curtime=1;

filename_pk2 = 'trash';
fout = fopen(get_output_file_name(filename_pk2), 'w');
[nliga_return] = nliga_returns( 30, geo, multi_mesh, mat, [], [], fout,eps0, k0 );

u = nliga_return.u;



% Generate the test-line in param-space
n_points = 100; % should be high enough to get a smooth translation of the Bishop Frame
ps_start = [0, 0.17];
ps_end = [0, 0.23]; % Expect transition around 0.2

points = [linspace(ps_start(1), ps_end(1), n_points)',  linspace(ps_start(2), ps_end(2), n_points)'];
s = linspace(0, 1, n_points);

p_param = [zeros(n_points, 1),  [linspace(0.5, 1, 0.5*n_points)'; linspace(0, 0.5, n_points*0.5)']];
p_param = [linspace(0, 1, 0.5*n_points)', ones(0.5 * n_points, 1); linspace(0, 1, 0.5*n_points)', zeros(0.5 * n_points, 1)];
%p_param = [zeros(n_points, 1), linspace(0, 1, n_points)'];
%p_param = [linspace(0, 1, n_points)', linspace(0, 1, n_points)'];
% Pre-allocate space to safe the line-results to

pk2_line = zeros(6,n_points);
x_line = zeros(3, n_points);
F_line = zeros(3,3,n_points);
vm_line = zeros(1, n_points);




% Generate element domain (scaled to actual corodinates)
mesh = multi_mesh;
for pi = 1:n_points
    % Identify element our point lies in / on
    p = points(pi, :);

    % Hard coded, but necessary
    %e = 36;
    %m = 3;
    if pi <= n_points/2%p(2) < 0.2
        % Element 80
        e = 80;
        m = 5;
        % Translate P into parameter space and compute Shape function
        % weights
    else
        e = 36;
        m = 3;
    end
            
    % Evaluate Stress and u at parameter point
    sctr = mesh.elNodeCnt(e,:);     % element control points index
    exyz = mesh.coords(sctr,:);  % element control points' coordinates
    nn = length(sctr);   % number of control points in the element
    nn3 = nn*3;          % degree of freedom of control points
    
    % now use globElNodeCnt instead
    %sctr = mesh.submeshes{1, m}.gloElNodeCnt(e, :);
    sctrB = zeros(1, nn3);      
    sctrB(1:3:nn3) = 3*sctr - 2;% displacement in x direction
    sctrB(2:3:nn3) = 3*sctr-1;  % displacement in y direction
    sctrB(3:3:end) = 3*sctr;    % displacement in z direction
    edsp = u(sctrB);
    edsp = reshape(edsp, 3, nn);


    % TODO: Remove after testing
    %[N,ders] = nurbs_derivatives( [p_param(pi, 1), p_param(pi, 2)],geo, mesh );
    [N, ders] = nurbs_derivatives( [p_param(pi, 1), p_param(pi, 2)],geo{1,m}, mesh.submeshes{1, m});
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

    % Store line data
    x_line(:, pi) = x;
    pk2_line(:, pi) = pk2;
    F_line(:, :, pi) = F;
    vm_line(:, pi) = von_mises(pk2cauchy(pk2, F)');
end 






figure
grid on
hold on
scatter(mesh.coords(:, 1), mesh.coords(:, 2), 20, "blue", "filled")
plot(x_line(1, :), x_line(2, :), "red")



% Scatter-plot the Von-Mises Stress over s
figure
grid on
hold on

yyaxis left
plot(s, vm_line, "blue")
scatter(s, vm_line, 20, "blue", "filled")
xlabel("Arc Length parameter")
ylabel("Von Mises")
title("Von Mises value along a vertical Line through Meshes 5 and 3")





% Visualize over arc-length

%mesh5 = build_iga_mesh(geo{1,5});

% Normales square as mesh
geo = geo_square([0,0], 0.4);
mesh5 = build_iga_mesh(geo);

k = mesh5.p+1;
r = 10;
r_params = [zeros(r, 1), linspace(0.7, 0.8, r)'];

Ns = zeros(r, k);
for ri = 1:r
    [N, ders] = nurbs_derivatives( [0, r_params(ri, 2)], geo, mesh5);
    Ns(ri, :) = N(1:k:k*k);
end

figure
grid on
hold on
xs = 1:k;
ax = gca;
for ri = 1:r
    currentColorIndex = ax.ColorOrderIndex;
    scatter(xs, Ns(ri, :), 20, 'filled', 'DisplayName', sprintf('%d', ri)); % 'filled' sieht oft schicker aus
    
    % 3. Index zurücksetzen
    ax.ColorOrderIndex = currentColorIndex;
    plot(xs, Ns(ri, :), 'HandleVisibility','off');
end
xlabel("Y_80")
ylabel("N Weight")
