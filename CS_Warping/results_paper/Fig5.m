% This script is used to test the equivalency of the PK2 Formulation
% against the PK1 Formulation. It includes:
%   - 3 Load Cases as possible test cases
%   - 2 cross-sectional shapes
%   - Different Material Indices
%   - Options to save and export results as a table
%   - Visualization options
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

% The following options may be chosen to reproduce the figures 5a and 5b:
% (The current preset corresponds to Fig 5a)
%           | cs_type   | eps0_case    | k0_case    | complex_case
%   Fig 5a  | "square"  | 0            | 0          | 1
%   Fig 5b  | "square"  | 0            | 0          | 1





clc;
clear;

%% Options

% Select the Loading case
eps0_case = 0;
k0_case = 0;
complex_case = 1;

% Select the Cross-section
crossectional_type = "square";
%crossectional_type = "circle";

% Select visualization options
vis_beam_effects = 1;
vis_beam_stiffnesses = 1;
vis_simple_beam_stiffnesses = 1;
show_errors = 1;

% Select if files should be generated
safe_files = 0;



% Select material modeling index -> Execute testcase for every index
%   - PK2 Formulations (110-119)
%   - PK1 Formulations (10-19)
index_SVK_pk1 = 14; % Saint-Venant Kirchhoff with PK1
index_SVK_pk2 = 114; % Saint-Venant Kirchhoff with PK2

index_MR_pk1 = 11; % Compressible Mooney-Rivling with PK1
index_MR_pk2 = 111; % Compressible Mooney-Rivling with PK2

index_NH_pk1 = 10;  % Compressible Neo-Hooke with PK1
index_NH_pk2 = 110;  % Compressible Neo-Hooke with PK2

indexes = [index_SVK_pk1, index_SVK_pk2, index_MR_pk1, index_MR_pk2, index_NH_pk1, index_NH_pk2];

% Program options
append_older_dataset = 0;


if eps0_case == 1
    % Z-Axial Strain prescribed
    dir = 3;
    eps0i = 0:0.01:0.1;
    l = length(eps0i);
    loadcase = zeros(6,l);
    loadcase(dir, :) = eps0i;

    save_tables_eps0= 1;
    save_tables_k03 = 0;
    x_axis_txt = 'Axial Strain $\epsilon_{3}$';
elseif k0_case == 1
    % Z-Axial Twist prescribed
    dir = 6;
    k0i = 0:0.01:0.1;
    l = length(k0i);
    loadcase = zeros(6,l);
    loadcase(dir, :) = koi;

    save_tables_eps0= 0;
    save_tables_k03 = 1;
    x_axis_txt = 'Axial Twist $\kappa_{3}$';
elseif complex_case == 1
    % Multiaxial loading case
    % Make 'substep' Increments of loading from Null-Load to Max-Load
    eps0_max = [0.02, 0.03, 0.1]';
    k0_max = [0.01,0.02,0.02]';
    
    n_substeps = 10; % Number of Substeps
    substep_range = 0:n_substeps;
    substep_range = substep_range / n_substeps; % Normalise Substep vector
    l = length(substep_range);
    loadcase = [eps0_max;k0_max] .* substep_range;
    
    x_axis_txt = 'Loading Increment factor $\lambda$';
end




% Build geometrical model
if crossectional_type == "square"
    cs_options = {};
    cs_options.RefinementX = 9;
    cs_options.RefinementY = 9;
    cs_options.show_plot = 0;
    plate = geo_square([0,0], 1, cs_options);
elseif crossectional_type == "circle"
    cs_options = {};
    cs_options.Refinement = 3;
    cs_options.show_plot = 0;
    plate = geo_circle_with_square( [0,0], 1, 0.6, cs_options);
else
    error("No matching crossectional shape selected")
end


% Enforce displacement boundary conditions 
dbc = [];        
tbc = [];
tol = 1e-8;
dof = 3;
eltype = 30;
titles = ["SVK (PK1)", "SVK (PK2)", ...
    "Mooney-Rivlin (PK1)", "Mooney-Rivlin (PK2)", ...
    "Neo-Hook (PK1)", "Neo-Hook (PK2)"];
colors = ["blue", "blue", ...
    "green", "green", ...
    "red", "red"];
shapes = ["square", "square", "^", "^", "o", "o"];

% Store beam forces and stiffnesses
l_start = 1;

k = length(indexes); % Number of Indizees to compute

% Build iga mesh structure
if crossectional_type == "circle"
    % execute custom multi-mesh assembly
    mesh = build_omesh(plate);
else 
    mesh = build_iga_mesh( plate );
end

% Preallocate space to store computed variables
if append_older_dataset
    % load older nl_data dataset
    load("PK_Comparison.mat");
    l_start = size(nl_data.n0, 3) + 1;
else
    nl_data.compute_time = zeros(l, k);
    nl_data.coords_def = zeros(3,mesh.nCpts,l, k);
    nl_data.n0 = zeros(3,1,l, k);
    nl_data.m0 = zeros(3,1,l, k);
    nl_data.C0 = zeros(6,6,l, k);
    s_u = 3*mesh.nCpts+6;
    nl_data.u = zeros(s_u, l, k);
    nl_data.k = zeros(s_u, s_u, l, k);
    nl_data.r = zeros(20, l, k); % Reserve 20 Rows for residuals for a load-step
end



% Compute the beam forces and stiffnesses for
%   Loading case
%   Material Model
%   PK1 / PK2 Formulation
for i = 1:l

    eps0 = loadcase(1:3, i);
    k0 = loadcase(4:6, i);

    l_index = l_start + (i - 1);
    disp("(" + num2str(i / l * 100) + " %) Computing for eps0/k0: " + num2str([eps0; k0]'))
    for j = 1:k
        
        % Retrieve material properties
        mat = default_mat();
        mat.index = indexes(j);
        filename = 'CSWP_PK_Comparison';
        fname = get_output_file_name(filename);
        fout = fopen(fname,'w'); 
        
        % Start Compute timer
        tic;
        nl_returns = nliga_returns( eltype, plate, mesh, mat, dbc, tbc, fout, eps0, k0);
        tc = toc;
        nl_data.compute_time(l_index, j) = tc;
        
        disp("\t (" + (i/l*100) + "%)... Done " + titles(j) + " (" + tc + ")");
          
        fclose(fout);

        % Store the data
        nl_data.coords_def(:, :, l_index, j) = nl_returns.coords_def';
        nl_data.n0(:, :, l_index, j) = nl_returns.n0;
        nl_data.m0(:, :, l_index, j) = nl_returns.m0;
        nl_data.C0(:, :, l_index, j) = nl_returns.C0;
        nl_data.u(:, l_index, j) = nl_returns.u;
        nl_data.k(:, :, l_index, j) = nl_returns.k;
    end
end

%% Export the data as csv with header

if eps0_case == 1 || k0_case == 1
    xx = loadcase(dir, :);
elseif complex_case == 1
    xx = substep_range;
end

% Data Extract for axial stretching
n3_SVK_pk1 = reshape(nl_data.n0(3,1,:,1), [], 1);
n3_SVK_pk2 = reshape(nl_data.n0(3,1,:,2), [], 1);
n3_MR_pk1 = reshape(nl_data.n0(3,1,:,3), [], 1);
n3_MR_pk2 = reshape(nl_data.n0(3,1,:,4), [], 1);
n3_NH_pk1 = reshape(nl_data.n0(3,1,:,5), [], 1);
n3_NH_pk2 = reshape(nl_data.n0(3,1,:,6), [], 1);

C33_SVK_pk1 = reshape(nl_data.C0(3,3,:,1), [], 1);
C33_SVK_pk2 = reshape(nl_data.C0(3,3,:,2), [], 1);
C33_MR_pk1 = reshape(nl_data.C0(3,3,:,3), [], 1);
C33_MR_pk2 = reshape(nl_data.C0(3,3,:,4), [], 1);
C33_NH_pk1 = reshape(nl_data.C0(3,3,:,5), [], 1);
C33_NH_pk2 = reshape(nl_data.C0(3,3,:,6), [], 1);

% Data extract for torsional flexing
m3_SVK_pk1 = reshape(nl_data.m0(3,1,:,1), [], 1);
m3_SVK_pk2 = reshape(nl_data.m0(3,1,:,2), [], 1);
m3_MR_pk1 = reshape(nl_data.m0(3,1,:,3), [], 1);
m3_MR_pk2 = reshape(nl_data.m0(3,1,:,4), [], 1);
m3_NH_pk1 = reshape(nl_data.m0(3,1,:,5), [], 1);
m3_NH_pk2 = reshape(nl_data.m0(3,1,:,6), [], 1);

C66_SVK_pk1 = reshape(nl_data.C0(6,6,:,1), [], 1);
C66_SVK_pk2 = reshape(nl_data.C0(6,6,:,2), [], 1);
C66_MR_pk1 = reshape(nl_data.C0(6,6,:,3), [], 1);
C66_MR_pk2 = reshape(nl_data.C0(6,6,:,4), [], 1);
C66_NH_pk1 = reshape(nl_data.C0(6,6,:,5), [], 1);
C66_NH_pk2 = reshape(nl_data.C0(6,6,:,6), [], 1);


% save data tables
if safe_files && eps0_case
    T = array2table([xx', ...
        n3_SVK_pk1, n3_SVK_pk2, ...
        n3_MR_pk1, n3_MR_pk2, ...
        n3_NH_pk1, n3_NH_pk2, ...
        C33_SVK_pk1, C33_SVK_pk2, ...
        C33_MR_pk1, C33_MR_pk2, ...
        C33_NH_pk1, C33_NH_pk2]);
    T.Properties.VariableNames(1:13) = {'v0_3', ...
        'n3_SVK_pk1','n3_SVK_pk2', ...
        'n3_MR_pk1', 'n3_MR_pk2', ...
        'n3_NH_pk1', 'n3_NH_pk2', ...
        'C33_SVK_pk1', 'C33_SVK_pk2', ...
        'C33_MR_pk1', 'C33_MR_pk2', ...
        'C33_NH_pk1', 'C33_NH_pk2'};
    writetable(T,'PAPER_v03_N3_C33.csv')
end

if safe_files && k0_case
    T = array2table([xx', ...
        m3_SVK_pk1, m3_SVK_pk2, ...
        m3_MR_pk1, m3_MR_pk2, ...
        m3_NH_pk1, m3_NH_pk2, ...
        C66_SVK_pk1, C66_SVK_pk2, ...
        C66_MR_pk1, C66_MR_pk2, ...
        C66_NH_pk1, C66_NH_pk2]);
    T.Properties.VariableNames(1:13) = {'k0_3', ...
        'm3_SVK_pk1','m3_SVK_pk2', ...
        'm3_MR_pk1', 'm3_MR_pk2', ...
        'm3_NH_pk1', 'm3_NH_pk2', ...
        'C66_SVK_pk1', 'C66_SVK_pk2', ...
        'C66_MR_pk1', 'C66_MR_pk2', ...
        'C66_NH_pk1', 'C66_NH_pk2'};
    writetable(T,'PAPER_k03_m3_C66.csv')
end

if safe_files && complex_case
    T = array2table([xx', ...
        n3_SVK_pk1, n3_SVK_pk2, ...
        n3_MR_pk1, n3_MR_pk2, ...
        n3_NH_pk1, n3_NH_pk2, ...
        m3_SVK_pk1, m3_SVK_pk2, ...
        m3_MR_pk1, m3_MR_pk2, ...
        m3_NH_pk1, m3_NH_pk2, ...
        C33_SVK_pk1, C33_SVK_pk2, ...
        C33_MR_pk1, C33_MR_pk2, ...
        C33_NH_pk1, C33_NH_pk2, ...
        C66_SVK_pk1, C66_SVK_pk2, ...
        C66_MR_pk1, C66_MR_pk2, ...
        C66_NH_pk1, C66_NH_pk2]);
    T.Properties.VariableNames(1:25) = {'lambda', ...
        'n3_SVK_pk1','n3_SVK_pk2', ...
        'n3_MR_pk1', 'n3_MR_pk2', ...
        'n3_NH_pk1', 'n3_NH_pk2', ...
        'm3_SVK_pk1','m3_SVK_pk2', ...
        'm3_MR_pk1', 'm3_MR_pk2', ...
        'm3_NH_pk1', 'm3_NH_pk2', ...
        'C33_SVK_pk1', 'C33_SVK_pk2', ...
        'C33_MR_pk1', 'C33_MR_pk2', ...
        'C33_NH_pk1', 'C33_NH_pk2', ...
        'C66_SVK_pk1', 'C66_SVK_pk2', ...
        'C66_MR_pk1', 'C66_MR_pk2', ...
        'C66_NH_pk1', 'C66_NH_pk2'};
    writetable(T,'PAPER_Complex_Loading.txt')
end

%% Visualisation of Results

if show_errors == 1
    % Compare the Errors in C0, N0 and M0 Computation between First and
    % second index (PK1 vs PK2)
    disp("Error C0 (PK1 - PK2): ")
    disp(nl_data.C0(:,:,1,1) - nl_data.C0(:,:,1,2));
    disp("Error N0 (PK1 - PK2): ")
    disp(nl_data.n0(:,:,1,1) - nl_data.n0(:,:,1,2));
    disp("Error M0: (PK1 - PK2): ")
    disp(nl_data.m0(:,:,1,1) - nl_data.m0(:,:,1,2));
end


if vis_beam_effects == 1
    % Normal Force
    beam_figure_axial=figure(1);
    set(beam_figure_axial,'name','Beam Forces Normal Z');
    grid on;
    hold on;
    title(["Normal Beam Forces over ", x_axis_txt], 'interpreter', 'latex')
    xlabel(x_axis_txt, 'interpreter', 'latex')
    ylabel("Normal Force in [kN]", 'interpreter', 'latex')
    for index = 1:k
        if mod(index, 2) == 1 % PK1
            yy_pk1 = reshape(nl_data.n0(3,1,:,index), 1, []);
            scatter(xx, yy_pk1, 40, "filled", colors(index), shapes(index), DisplayName=titles(index))
            plot(xx, yy_pk1, colors(index), LineStyle="--", HandleVisibility="off")
        elseif mod(index, 2) == 0 % PK2
            yy_pk2 = reshape(nl_data.n0(3,1,:,index),1, []);
            scatter(xx, yy_pk2, 100, colors(index), shapes(index), DisplayName=titles(index))
        end
        legend() 
    end
    
    % Normalised Normal Force
    beam_figure_axial=figure(5);
    set(beam_figure_axial,'name','Normalised Beam Forces Normal Z');
    grid on;
    hold on;
    title(["Normalised Normal Beam Forces over ", x_axis_txt], 'interpreter', 'latex')
    xlabel(x_axis_txt, 'interpreter', 'latex')
    ylabel("Normalised Normal Force in [%]", 'interpreter', 'latex')
    for index = 1:k
        if mod(index, 2) == 1 % PK1
            yy_pk1 = reshape(nl_data.n0(3,1,:,index), 1, []);
            scatter(xx, yy_pk1 ./ max(yy_pk1), 40, "filled", colors(index), shapes(index), DisplayName=titles(index))
            plot(xx, yy_pk1 ./ max(yy_pk1), colors(index), LineStyle="--", HandleVisibility="off")
        elseif mod(index, 2) == 0 % PK2
            yy_pk2 = reshape(nl_data.n0(3,1,:,index),1, []);
            scatter(xx, yy_pk2 ./ max(yy_pk2), 100, colors(index), shapes(index), DisplayName=titles(index))
        end
        legend('Location', 'northwest') 
    end
    
    % Shear Force
    beam_figure_shear=figure(2);
    set(beam_figure_shear,'name','Beam Forces Shear X');
    grid on;
    hold on;
    title(["Shear Beam Forces over ", x_axis_txt], 'interpreter', 'latex')
    xlabel(x_axis_txt, 'interpreter', 'latex')
    ylabel("Shear Force in [kN]", 'interpreter', 'latex')
    for index = 1:k
        if mod(index, 2) == 1 % PK1
            yy_pk1 = reshape(nl_data.n0(1,1,:,index), 1, []);
            scatter(xx, yy_pk1, 40, "filled", colors(index), shapes(index), DisplayName=titles(index))
            plot(xx, yy_pk1, colors(index), LineStyle="--", HandleVisibility="off")
        elseif mod(index, 2) == 0 % PK2
            yy_pk2 = reshape(nl_data.n0(1,1,:,index),1, []);
            scatter(xx, yy_pk2, 100, colors(index), shapes(index), DisplayName=titles(index))
        end
        legend('Location', 'northwest') 
    end
    
    % Torsion Moment
    beam_figure_torsion=figure(3);
    set(beam_figure_torsion,'name','Beam Moments Torsion Z');
    grid on;
    hold on;
    title(["Torsional Beam Moment over ", x_axis_txt], 'interpreter', 'latex')
    xlabel(x_axis_txt, 'interpreter', 'latex')
    ylabel("Beam Moment in [Nm]", 'interpreter', 'latex')
    for index = 1:k
        if mod(index, 2) == 1 % PK1
            yy_pk1 = reshape(nl_data.m0(3,1,:,index), 1, []);
            scatter(xx, yy_pk1, 40, "filled", colors(index), shapes(index), DisplayName=titles(index))
            plot(xx, yy_pk1, colors(index), LineStyle="--", HandleVisibility="off")
        elseif mod(index, 2) == 0 % PK2
            yy_pk2 = reshape(nl_data.m0(3,1,:,index),1, []);
            scatter(xx, yy_pk2, 100, colors(index), shapes(index), DisplayName=titles(index))
        end
        legend('Location', 'northwest')  
    end
    
    % Bending Moment
    beam_figure_bending=figure(4);
    set(beam_figure_bending,'name','Beam Moments Bending X');
    grid on;
    hold on;
    title(["Bending Beam Moment over ", x_axis_txt], 'interpreter', 'latex')
    xlabel(x_axis_txt, 'interpreter', 'latex')
    ylabel("Beam Moment in [Nm]", 'interpreter', 'latex')
    for index = 1:k
        if mod(index, 2) == 1 % PK1
            yy_pk1 = reshape(nl_data.m0(1,1,:,index), 1, []);
            scatter(xx, yy_pk1, 40, "filled", colors(index), shapes(index), DisplayName=titles(index))
            plot(xx, yy_pk1, colors(index), LineStyle="--", HandleVisibility="off")
        elseif mod(index, 2) == 0 % PK2
            yy_pk2 = reshape(nl_data.m0(1,1,:,index),1, []);
            scatter(xx, yy_pk2, 100, colors(index), shapes(index), DisplayName=titles(index))
        end
        legend('Location', 'northwest') 
    end
end

if vis_simple_beam_stiffnesses == 1
    % Beam Stiffness
    
    % Normal Stiffness C0(3,3)
    beam_figure_axial=figure(6);
    set(beam_figure_axial,'name','Beam Stiffness Normal Z');
    grid on;
    hold on;
    title(["Normal Beam Stiffness over ", x_axis_txt], 'interpreter', 'latex')
    xlabel(x_axis_txt, 'interpreter', 'latex')
    ylabel("Axial Beam Stiffness in [kN]", 'interpreter', 'latex')
    for index = 1:k
        if mod(index, 2) == 1 % PK1
            yy_pk1 = reshape(nl_data.C0(3,3,:,index), 1, []);
            scatter(xx, yy_pk1, 40, "filled", colors(index), shapes(index), DisplayName=titles(index))
            plot(xx, yy_pk1, colors(index), LineStyle="--", HandleVisibility="off")
        elseif mod(index, 2) == 0 % PK2
            yy_pk2 = reshape(nl_data.C0(3,3,:,index),1, []);
            scatter(xx, yy_pk2, 100, colors(index), shapes(index), DisplayName=titles(index))
        end
        legend() 
    end
end



if vis_beam_stiffnesses == 1
    % Multiplot for all stiffness entries
    figupperleft = figure('Name', 'C0 Entries upper left quadrant', 'NumberTitle', 'off');
    figupperright = figure('Name', 'C0 Entries upper right quadrant', 'NumberTitle', 'off');
    figlowerleft = figure('Name', 'C0 Entries lower left quadrant', 'NumberTitle', 'off');
    figlowerright = figure('Name', 'C0 Entries lower right quadrant', 'NumberTitle', 'off');
    
    figs = [figupperleft, figupperright, figlowerleft, figlowerright];
    fig_titles = ["upper left", "upper right", "lower left", "lower right"];
    offsets = [0 0; 0 3; 3 0; 3 3];


    for fig_index = 1:4
        fig = figs(fig_index);
        t = tiledlayout(fig, 3,3);
        title(t, "Beam Stiffness Entry (" + fig_titles(fig_index) + " quadrant) over " + x_axis_txt, 'interpreter', 'latex')
        xlabel(t, x_axis_txt, 'interpreter', 'latex')
        ylabel(t, "Beam Stiffness in [kN]", 'interpreter', 'latex')

        for row = 1:3
            for col = 1:3
                nexttile(t)
                row_col_data = offsets(fig_index, :) + [row col];
                grid on
                hold on
                for index = 1:k
                    
                    if mod(index, 2) == 1 % PK1
                        yy_pk1 = reshape(nl_data.C0(row_col_data(1),row_col_data(2),:,index),1, []);
                        scatter(xx, yy_pk1, 40, "filled", colors(index), shapes(index), DisplayName=titles(index))
                        plot(xx, yy_pk1, colors(index), LineStyle="--", HandleVisibility="off")
                    elseif mod(index, 2) == 0 % PK2
                        yy_pk2 = reshape(nl_data.C0(row_col_data(1),row_col_data(2),:,index),1, []);
                        scatter(xx, yy_pk2, 100, colors(index), shapes(index), DisplayName=titles(index))
                    end
                    legend() 
                end
            end
        end
    end
end
