% This script is used to validate the PK1 and PK2 Formulations against
% results by Arora (see (1) for more details).
% ------------------------------------------------------------------------ 
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

% setup the loading case
eps0i = 0;
k0i = 0:0.05:0.5;
n_datapoints = max(length(eps0i), length(k0i));
d = 3;
append_older_dataset = 0;

eps0s = zeros(3, n_datapoints);
eps0s(d, :) = eps0s(d, :) + eps0i;


k0s = zeros(3, n_datapoints);
k0s(3, :) = k0i;
crossectional_type = "square";
x_axis_txt = "$k_{0_{" + num2str(d) + "}}$\,(Axial Twist)";

% Options
export_to_csv = 1;
visualize_validation = 1;

% Build geometrical model (See Arora)
plate =  geo_rectangle( [0,0], 1, 2);


% Enforce displacement boundary conditions 
dbc = [];        % dbc = [node index, node dof, prescribed displacement]
tbc = [];
tol = 1e-8;

% Determine material properties
% Note that definition 'mat' is different for nonlinear materials, the first number
% is always the index to define material categories
% index - [20-40) - elastoplastic material

dof = 3;
eltype = 30;
index_SVK_pk1 = 14; % Saint-Vernant Kirchhoff with PK1
index_SVK_pk2 = 114; % Saint-Vernant Kirchhoff with PK2

index_MR_pk1 = 11; % Compressible Mooney-Rivling with PK1
index_MR_pk2 = 111; % Compressible Mooney-Rivling with PK2

index_NH_pk1 = 10;  % Compressible Neo-Hooke with PK1
index_NH_pk2 = 110;  % Compressible Neo-Hooke with PK2

%indexes = [index_SVK_pk1, index_SVK_pk2, index_MR_pk1, index_MR_pk2, index_NH_pk1, index_NH_pk2];
indexes = [index_SVK_pk1, index_SVK_pk2];
titles = ["SVK (PK1)", "SVK (PK2)", ...
    "Mooney-Rivlin (PK1)", "Mooney-Rivlin (PK2)", ...
    "Neo-Hook (PK1)", "Neo-Hook (PK2)"];
colors = ["blue", "blue", ...
    "green", "green", ...
    "red", "red"];
shapes = ["square", "square", "^", "^", "o", "o"];

% Store beam forces and stiffnesses
l = size(eps0s, 2); % Number of V0(3) positions to compute
l_start = 1;

k = length(indexes); % Number of PK Formulations to compute

if append_older_dataset
    % load older nl_data dataset
    load("PK_Comparison.mat");
    l_start = size(nl_data.n0, 3) + 1;
else
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
    eps0_i = eps0s(d, i);
    k0_i = k0s(d, i);

    % expect onset of "flow" at 1.95 (???)
    eps0 = eps0s(:, i);
    k0 = k0s(:, i);
    
    % Build iga mesh structure
    mesh = build_iga_mesh( plate );
    
    % Build four edges
    curve = extract_iga_boundary(mesh);
    
    l_index = l_start + (i - 1);
    nl_data.eps0 (:,:, l_index) = eps0;
    disp("(" + num2str(i / l * 100) + " %) Computing for eps0/k0: " + num2str([eps0; k0]'))
    for j = 1:k
        
        % Retrieve material properties
        % Changed Material properties to fit with Arora Diss

        mat = Arora_mat();
        mat.index = indexes(j);

        filename = 'ANALYSIS_Arora_Validation';
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


% Export generated Data to .csv files
if export_to_csv == 1
    xx = reshape(nl_data.k0(d,:),1,[]);

    % Data extract for torsional Moment
    C66_SVK_pk1 = reshape(nl_data.C0(6,6,:,1), [], 1);
    C66_SVK_pk2 = reshape(nl_data.C0(6,6,:,2), [], 1);

    % Normalise data
    C66_SVK_pk1 = C66_SVK_pk1 ./ C66_SVK_pk1(1);
    C66_SVK_pk2 = C66_SVK_pk2 ./ C66_SVK_pk2(1);

    % save data table
    T = array2table([xx', ...
        C66_SVK_pk1, C66_SVK_pk2]);
    T.Properties.VariableNames(1:3) = {'k0_3', ...
        'c66_SVK_PK1','c66_SVK_PK2'};
    writetable(T,'PAPER_Validation_Arora_Us.csv')
end



if visualize_validation == 1
    % Load in Arora-Data
    Arora_k = readtable('digitized_datasets/AroraStiffness.txt');
    
    % "Normalised" Torsional Stiffness
    beam_figure_axial=figure(105);
    set(beam_figure_axial,'name','Validation Arora');
    grid on;
    hold on;
    ylabel("Normalised Torsional Stiffness")
    xx = reshape(nl_data.k0(3,:),1,[]);
    xlabel("K03")
    for index = 1:k
        if mod(index, 2) == 1 % PK1
            yy_pk1 = reshape(nl_data.C0(6,6,:,index), 1, []);
            scatter(xx, yy_pk1 ./ yy_pk1(1), 40, "filled", colors(index), shapes(index), DisplayName=titles(index))
            plot(xx, yy_pk1 ./ yy_pk1(1), colors(index), LineStyle="--", HandleVisibility="off")
        elseif mod(index, 2) == 0 % PK2
            yy_pk2 = reshape(nl_data.C0(6,6,:,index),1, []);
            scatter(xx, yy_pk2 ./ yy_pk2(1), 100, colors(index), shapes(index), DisplayName=titles(index))
        end
    
        legend() 
    end
    
    % Plot Arora Results
    scatter(Arora_k.k', Arora_k.SVK', "red", DisplayName="Arora")
end