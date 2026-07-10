%% FOR TESTING AND DEMONSTRATION PURPOSES
% This Demonstration Environment showcases the computation of the beam 
% forces and beam stiffnesses
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


% Define the loadcase (10% axial stretch)
% Attention: Always use upright vectors
eps0 = [0, 0, 0.1]';
k0 = [0, 0, 0]';


eps0 = [0.01, 0.04, 0.06]';
k0 = [-0.05, 0.04, 0.09]';


% Define Element type, Safe File and Boundary conditions
eltype = 30; % 30- CSWP element
filename = 'TEST_newStress';
fname = get_output_file_name(filename);
fout = fopen(fname,'w'); 
dbc =[]; % Dirichlet boundary conditions
tbc=[]; % Von Neumann boundary conditions

% Define Material, Material Model Type, Geometry, Mesh
mat = default_mat();
mat.index = 114; % SVK with PK2 / Alternatives see "default_mat()"

use_square_mesh = 1;

if use_square_mesh
    % Use Square Mesh
    geo = geo_square([0,0], 1);
    mesh = build_iga_mesh( geo );
else
    % Use the circular O-Mesh Formulation
    geo = geo_circle_with_square([0,0], 1, 0.4);
    mesh = build_omesh(plate_circle);
end


if use_scaled_fullaxial_loading:
    num_loading_steps = 10;
    loadings = linspace(0, 1, num_loading_steps) .* [eps0; k0];
else
    % All possible Uniaxial Loadings
end



n0s = zeros(3, num_loading_steps);
m0s = zeros(3, num_loading_steps);
C0s = zeros(6,6,num_loading_steps);
time_old = zeros(num_loading_steps);


n0s_new = zeros(3, num_loading_steps);
m0s_new = zeros(3, num_loading_steps);
C0s_new = zeros(6,6,num_loading_steps);
time_new = zeros(num_loading_steps);

for loading_scenario = 1:num_loading_steps
    
    % Extract loading
    eps0 = loadings(1:3, loading_scenario);
    k0 = loadings(4:6, loading_scenario);

    % Solve the NLIGA simulation for the Square
    %nl_return = nliga_returns(eltype, geo, mesh, mat, dbc, tbc, fout, eps0, k0);
    
    % Solve the NLIGA simullation for the O-mesh
    nl_return = nliga_returns(eltype, geo, mesh, mat, dbc, tbc, fout, eps0, k0);
    u = nl_return.u; % Solution displacement vector
    k = nl_return.k; % Solution stiffness matrix
    
    % Extract old stresses and stiffness
    n0s(:, loading_scenario) = nl_return.n0;
    m0s(:, loading_scenario) = nl_return.m0;
    C0s(:,:,loading_scenario) = nl_return.C0;
    time_old(loading_scenario) = nl_return.nmc_time;

    % Extract new stresses and stiffness
    n0s_new(:, loading_scenario) = nl_return.n0_new;
    m0s_new(:, loading_scenario) = nl_return.m0_new;
    C0s_new(:,:,loading_scenario) = nl_return.C0_new;
    time_new(loading_scenario) = nl_return.nmc_time_new;
end

% Visualize the absolute and relative errors over load scenarios
xx = 1:num_loading_steps;

y_n_abs = sqrt(sum((n0s - n0s_new).^2, 1));
y_m_abs = sqrt(sum((m0s - m0s_new).^2, 1));
y_c0_abs = sqrt(sum(sum((C0s - C0s_new).^2, 1), 2));

figure('Position', [100, 100, 1000, 800])
sgtitle("Square Mesh | Linear scaling of full-axial loading case")


subplot(2, 2, 1)
plot(xx, y_n_abs, 'r-', 'LineWidth', 1.5)
hold on
scatter(xx, y_n_abs, 20, 'r', 'filled')
grid on
xlabel('Loading Step')
ylabel('Absolute Error (L2 Norm)')
title('Force Resultant Error (n0)')

subplot(2, 2, 2)
plot(xx, y_m_abs, 'b-', 'LineWidth', 1.5)
hold on
scatter(xx, y_m_abs, 20, 'b', 'filled')
grid on
xlabel('Loading Step')
ylabel('Absolute Error (L2 Norm)')
title('Moment Resultant Error (m0)')

subplot(2, 2, 3)
plot(xx, squeeze(y_c0_abs), 'g-', 'LineWidth', 1.5)
hold on
scatter(xx, squeeze(y_c0_abs), 20, 'g', 'filled')
grid on
xlabel('Loading Step')
ylabel('Absolute Error (Frobenius Norm)')
title('Stiffness Matrix Error (C0)')

subplot(2, 2, 4)
plot(xx, time_old(:, 1), 'k--', 'LineWidth', 1.5)
hold on
plot(xx, time_new(:, 1), 'm-', 'LineWidth', 1.5)
scatter(xx, time_old(:, 1), 20, 'k', 'filled')
scatter(xx, time_new(:, 1), 20, 'm', 'filled')
grid on
xlabel('Loading Step')
ylabel('Computation Time (s)')
title('Computation Time Comparison')
legend('Old Method', 'New Method', 'Location', 'best')