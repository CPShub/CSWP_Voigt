%% FOR TESTING AND DEMONSTRATION PURPOSES
% This Demonstration Environment shows the different visualization tools
% and options developed
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
clc;
clear;

recompute = 1;
filename1 = 'CSWP_visual_test1';
filename2 = 'CSWP_visual_test2';

fname1 = get_output_file_name(filename1);
fname2 = get_output_file_name(filename2);

if recompute
    % Set of Strain prescriptors for the first file
    v01 = [0, 0, 0]';
    k01 = [0.1, 0, 0.2]';

    % Set of Strain prescriptors for the second file
    v02 = [0.15, 0.15, 0]';
    k02 = [0, 0, 0]';
    
    
    % Define crossection and material
    cs_size = 1; % Unit Square
    square = geo_square([0, 0], cs_size, 0);
    mesh = build_iga_mesh(square);
    mat = default_mat();
    mat.index = 114; %SVK using PK2 Formulation

    
    % Constants for the nliga call
    eltype = 30;    % element type: 10 - plane strain element, 20 - solid element, 30 - CSWP element
    dbc =[];
    tbc=[];
    
    % Compute the first CSWP 
    fout1 = fopen(fname1,'w'); 
    disp("... Computing the first CSWP for eps0: [" + num2str(v01') + "] | k0: [" + num2str(k01') + "]")
    nliga_return1 = nliga_returns( eltype, square, mesh, mat, dbc, tbc, fout1 ,v01,k01);

    % Compute the second CSWP 
    fout2 = fopen(fname2,'w'); 
    disp("... Computing the second CSWP for eps0: [" + num2str(v02') + "] | k0: [" + num2str(k02') + "]")
    nliga_return2 = nliga_returns( eltype, square, mesh, mat, dbc, tbc, fout2 ,v02,k02);
end

%% Visualizing the CSWP-Results
% Visualize as a 2D flat color plot (version used in the Paper)
flag_flatcolor = 3; % Visualize the 3rd-component of displacement u
options.show_ticks = 1;
options.show_title = 1;
options.show_coords = 1;
options.show_cb_title = 1;
options.given_title = "Look, a title!";
plot_color_flat(flag_flatcolor, fname1, options);


% Visualize a direct comparison of the two loading cases using two
% flat-color visualizations
flag_comparison = 3; % Visualize the 3rd-component of displacement u
cs_coords_center = [-1, -1];
cs_A_text_pos = [0.4, 0.4];
cs_B_text_pos = [-0.6, -0.5];
options = {};
options.show_ticks = 0;
options.show_title = 1;
options.show_coords = {};
options.show_coords.flag = 1;
options.show_coords.center = cs_coords_center;
options.given_title = "Look, a title!";
options.A_text.pos = cs_A_text_pos;
options.A_text.text = "File 1";
options.A_text.font = 28;
options.A_text.color = "black";
options.B_text.pos = cs_B_text_pos;
options.B_text.text = "File 2";
options.B_text.font = 28;
options.B_text.color = "black";
options.loading_case = "u3";
options.cb_decimals = 3;

plot_color_flat_combined(flag_comparison, fname1, fname2, options);
