%% FOR TESTING AND DEMONSTRATION PURPOSES
% This Demonstration Environment showcases the computation of the 
% crossectional sensitivity to the loadcase variables [eps0, k0] by
% comparing it with a numerical approximation.
% ----------------------------------------
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




% Define the loadcase (10% axial stretch)
% Attention: Always use upright vectors

eps0 = 0; % Starting / Convergence Point
k0 = [0,0,0]';
y_t_integrals = zeros(1, 1);

% Define Options for visualization
plot_quivers = 1;
check_convergence = 1;

if check_convergence
    deps0 = [0, 0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001];
    eps0s = eps0 + deps0;
    steps = length(eps0s);
else
    deps = 0.05;
    steps = 6;
    eps0s = [eps0:deps:eps0+(steps-1)*deps];
end

% Define Element type, Safe File and Boundary conditions
eltype = 30; % 30- CSWP element
filename = 'DEMO_Beam_Effects';
fname = get_output_file_name(filename);
fout = fopen(fname,'w'); 
dbc =[]; % Dirichlet boundary conditions
tbc=[]; % Von Neumann boundary conditions

% Define Material, Material Model Type, Geometry, Mesh
mat = default_mat();
mat.index = 114; % SVK with PK2 / Alternatives see "default_mat()"
geo = geo_square([0,0], 1);
mesh = build_iga_mesh( geo );

% Allocate space for the positions and sensitivities over the time steps
xs = zeros(64, 3, steps);
ys = zeros(64, 3, steps);

for i = 1:steps
    % Compute K and u for small increments of eps0
    eps0 = [0 0 eps0s(i)]';
    
    % Reset the Geometries
    square = geo_square([0, 0], 1, 0);
    mesh = build_iga_mesh(square);
    
    % Solve the CSWP
    disp("... (" + num2str(i/steps * 100) + "%) Computing C0 for eps0: [0 0 " + num2str(eps0s(i)) + "]");
    nliga_return = nliga_returns( eltype, square, mesh, mat, dbc, tbc, fout ,eps0,k0);
    K = nliga_return.k;
    u = nliga_return.u;
    
    % Compute the crossectional sensitivity
    y_t = crossectional_stretch_sensitivity(square, mesh, mat, eps0, k0, u, K);
    
    % Compute single-value integral of y_t over area
    arrowdata = reshape(y_t(:, 3),3,[])';
    y_t_integrated = sum(arrowdata(1, :));

    % Store the data
    xs(:, :, i) = mesh.coords(:, [1,2,3]) + reshape(u(1:end-6), 3, [])';
    ys(:, :, i) = arrowdata;

    % Plot the crossectional stretch sensitivity
    if plot_quivers && i == 1
        fff = 1;
        poss = mesh.coords(:, [1,2,3]);
        figure
        hold on
        grid on
        scatter(poss(:, 1), poss(:, 2), 20, "red", "filled", 'DisplayName', 'Integration Point');
        quiver(poss(:, 1), poss(:, 2), fff * arrowdata(:, 1), fff * arrowdata(:, 2), 'DisplayName', 'Def. Senstivity');
        title("Axial Stretch: " + num2str(eps0s(i)));
        xlabel("X Axis")
        ylabel("Y Axis")
        legend()
        hold off
    end
end


% Approximate the sensitivity of displacement to the input variables (y)
% via the numerical difference dx / dv
if check_convergence
    % compute and plot the convergence of the error with decreasing dv
    dxs = xs(:,:,2:end) - xs(:,:,1);
    y_comp = zeros(size(dxs));
    for id = 1:length(deps0)-1
        y_comp(:,:,id) = dxs(:,:,id) ./ deps0(id+1);
    end

    y_err = ys(:,:,1:end-1) - y_comp;
    y_err_rel = y_err ./ y_comp;
    
    % plot the relative error over eps0 (Should converge to zero)
    mean_err = reshape(mean(y_err_rel(:,[1,2],:), [1, 2]), 1, []);
    
    figure
    loglog(deps0(2:end), abs(mean_err), '-ro');
    hold on
    grid on
    set(gca, 'XScale', 'log', 'YScale', 'log'); 
    set(gca, 'XDir', 'reverse');              
    xlabel('Axial Stretch stepsize (dv)')
    ylabel('|Rel. Error|')
    title('Convergence of relative error between analytical and numerical solution')
    xlim([min(deps0(2:end)) max(deps0(2:end))])
end