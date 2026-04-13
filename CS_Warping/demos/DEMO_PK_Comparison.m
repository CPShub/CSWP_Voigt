%% FOR TESTING AND DEMONSTRATION PURPOSES
% This Demonstration Environment showcases the equivalency of both the PK1
% and PK2 formulation of the CSWP and how both can be employed using the
% implementation
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
% notation", [Journal Name], [Year]. DOI: [DOI] 
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




clc;
clear;


% Define Element type, Safe File and Boundary conditions
eltype = 30; % 30- CSWP element
filename = 'DEMO_PK_Comparison';
fname = get_output_file_name(filename);
fout = fopen(fname,'w'); 
dbc =[]; % Dirichlet boundary conditions
tbc=[]; % Von Neumann boundary conditions

% Define Geometry, Mesh
geo = geo_square([0,0], 1);
mesh = build_iga_mesh( geo );
crossectional_type = "square";


% Set Options

% - Visualisation options
vis_beam_effects = 1;
vis_beam_stiffnesses = 1;
show_errors = 1;

% - Programm options
append_older_dataset = 0;
full_verification_check = 0;


if full_verification_check == 1
    % Executed load cases:
    %   - (1-3) Uniaxial Loading (eps1 || 2 || 3)
    %   - (4-6) Uniaxial Torsion (k01 || 2 || 3)
    %   - (7) Combined Axial Loading ( eps1 && 2 && 3)
    %   - (8) Combined Torsion (k01 && 2 && 3)
    %   - (9) Random Full loading (k01 && 2 && 3 && eps1 && 2 && 3)

    eps0s = zeros(3, 7);
    eps0s(3,:) = 1;
    k0s =zeros(3, 7);

    randset_eps = [0.04, 0.06, 0.05];
    randset_k0 = [0.08, 0.02, 0.1];

    eps0s(1, [1, 7]) = randset_eps(1);
    eps0s(2, [2, 7]) = randset_eps(2);
    eps0s(3, [3, 7]) = randset_eps(3);

    k0s(1, [4, 7]) = randset_k0(1);
    k0s(2, [5, 7]) = randset_k0(2);
    k0s(3, [6, 7]) = randset_k0(3);
else

    % You can play around with these to test different loading scenarios
    eps0i = 0:0.01:0.1;
    k0i = 0;
    n_datapoints = max(length(eps0i), length(k0i));
    d = 3;
    
    eps0s = zeros(3, n_datapoints);
    eps0s(d, :) = eps0s(d, :) + eps0i;
    k0s = zeros(3, n_datapoints);
    x_axis_txt = "$eps_{0_{3}}$\,(Axial Stretch)";
end



% Build geometrical model
if crossectional_type == "square"
    plate =  geo_square( [0,0], 1, 0);
elseif crossectional_type == "circle"
    plate = geo_circle( [0, 0], 1);
else
    error("No matching crossectional shape selected")
end

index_SVK_pk1 = 14; % Saint-Vernant Kirchhoff with PK1
index_SVK_pk2 = 114; % Saint-Vernant Kirchhoff with PK2

index_MR_pk1 = 11; % Compressible Mooney-Rivling with PK1
index_MR_pk2 = 111; % Compressible Mooney-Rivling with PK2

index_NH_pk1 = 10;  % Compressible Neo-Hooke with PK1
index_NH_pk2 = 110;  % Compressible Neo-Hooke with PK2

indexes = [index_SVK_pk1, index_SVK_pk2, index_MR_pk1, index_MR_pk2, index_NH_pk1, index_NH_pk2];
titles = ["SVK (PK1)", "SVK (PK2)", ...
    "Mooney-Rivlin (PK1)", "Mooney-Rivlin (PK2)", ...
    "Neo-Hook (PK1)", "Neo-Hook (PK2)"];
colors = ["blue", "blue", ...
    "green", "green", ...
    "red", "red"];
shapes = ["square", "square", "^", "^", "o", "o"];

% Store beam forces and stiffnesses
l = size(eps0s, 2); % Number of eps0(3) positions to compute
l_start = 1;

k = length(indexes); % Number of PK Formulations to compute

if append_older_dataset
    % load older nl_data dataset
    % Adapt this if necessary
    load("PK_Comparison.mat");
    l_start = size(nl_data.n0, 3) + 1;
else
    % Preallocate space for all relevant datasets
    nl_data.compute_time = zeros(l, k);
    nl_data.k0 = k0s;
    nl_data.eps0 = zeros(3,1,l);
    nl_data.coords_def = zeros(3,64,l, k);
    nl_data.n0 = zeros(3,1,l, k);
    nl_data.m0 = zeros(3,1,l, k);
    nl_data.C0 = zeros(6,6,l, k);
    nl_data.u = zeros(198, l, k);
    nl_data.k = zeros(198, 198, l, k);
end

% Compute the beam forces and stiffnesses for
%   Loading case
%   Material Model
%   PK1 / PK2 Formulation
for i = 1:l
    eps0 = eps0s(:, i);
    k0 = k0s(:, i);
    
    % Build iga mesh structure
    mesh = build_iga_mesh( plate );
    
    l_index = l_start + (i - 1);
    nl_data.eps0 (:,:, l_index) = eps0;
    disp("(" + num2str(i / l * 100) + " %) Computing for eps0/k0: " + num2str([eps0; k0]'))
    for j = 1:k
        
        % Retrieve material properties and set material index to 
        % indicate the PK Formulation and Material model to be used
        
        % This may be changed, see material models in functions/materials/
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

%% Visualisation of Results

if show_errors == 1
    disp("Error C0: ")
    disp(nl_data.C0(:,:,1,1) - nl_data.C0(:,:,1,2));
    disp("Error N0: ")
    disp(nl_data.n0(:,:,1,1) - nl_data.n0(:,:,1,2));
    disp("Error M0: ")
    disp(nl_data.m0(:,:,1,1) - nl_data.m0(:,:,1,2));
end


if full_verification_check
    % Only look at the loading index
    xx = 1:l;
else
    % Only look at the single changed loading variable
    xx = reshape(nl_data.eps0(d,1,:),1,[]);
end

% Compare the beam forces and moments resulting from both formulations
if vis_beam_effects == 1
    % Beam Forces
    beam_figure_axial=figure(1);
    set(beam_figure_axial,'name','Beam Forces Normal Z');
    grid on;
    hold on;
    title(["Normal Beam Forces over ", x_axis_txt], 'interpreter', 'latex')
    xlabel(x_axis_txt, 'interpreter', 'latex')
    ylabel("Normal Force N in Newton")
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
    ylabel("Normalised Normal Force")
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
    ylabel("Shear Force N in Newton")
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
    ylabel("Beam Moment in Nm")
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
    ylabel("Beam Moment in Nm")
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
        ylabel(t, "Beam Stiffness in $\frac{N}{mm^2}$", 'interpreter', 'latex')

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