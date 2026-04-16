% This script is used to generate combined gradient-colored plots of
% different displacement and stress derivatives over the cross-sections of
% two (different) output files, thus allowing an easy and direct visual
% comparison.
% Representable Options:
% Cross-Section
%   - Circular Cross-Section
%   - Square Cross-Section
% Dataset
%   - U3 Displacement Distribution
%   - Von-Mises Stress Distribution
% Loading Case
%   - Simple X-Axis Bending
%   - Multi-Axial (Full) Loading
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

% The following options may be chosen to reproduce the figures 4a and 4b:
% (The current preset corresponds to Fig 4a)
%           | cs_type   | loading_case      | display_type
%   Fig 4a  | "square"  | "full"            | "u3"
%   Fig 4b  | "circle"  | "full"            | "u3"

%   Fig 6a  | "square"  | "simple"          | "u3"
%   Fig 6b  | "circle"  | "simple"          | "u3"
%   Fig 6c  | "square"  | "simple"          | "vm"
%   Fig 6d  | "circle"  | "simple"          | "vm"



% To orient the dividing line diagonally, set the following value to 0
use_vertical_divide = 1;


% Select cross-section, loading case and visualized data
%cs_type = "circle";
cs_type = "square";
%cs_type = "circle_square";

loading_case = "full"; % Multi-Axial Loading case
%loading_case = "simple";% Uni-Axial X-Shear


display_type = "u3";    % u3 displacement component
%display_type = "vm";   % von-Mises Stress

recompute_files = 0;    % display_type may be changes without recomputing
save_file = 0;  


%%
options = {};

if use_vertical_divide
    options.diagonal.A = [0,2];
    options.diagonal.B = [0,-2];
    cs_A_text_pos = [0.43, -0.02];
    cs_B_text_pos = [-0.77, -0.02];
else
    options.diagonal.A = [-2,2];
    options.diagonal.B = [2,-2];
    cs_A_text_pos = [0.4, 0.4];
    cs_B_text_pos = [-0.6, -0.5];
end


cs_coords_center = [-1, -1];

% Cross-Section
if cs_type == "circle"
    plate = geo_circle( [0, 0], 1);
    savefile_cs = '_Circle.jpg';
elseif cs_type == "square"
    plate = geo_square( [0,0], 1, 0);
    savefile_cs = '_Square.jpg';
elseif cs_type == "circle_square"
    plate = geo_circle_with_square( [0,0], 1, 0.6);
    savefile_cs = "_CircleSquare.jpg";
end

% Dataset
if display_type == "u3"
    display_flag_text = '$U_{3}$ Displacement';
    display_flag = 3; % U3-Component
    savefile_display_type = '_U3';
else
    display_flag = 8; % Von Mises
    display_flag_text = 'Von-Mises Stress';
    savefile_display_type = '_VM';
end

% Loading Case
if loading_case == "full"
    % Full Loading
    eps0 = [0.02, 0.03, 0.1]';
    k0 = [0.01,0.02,0.02]';
    title_ = ['Multiaxial Loading - ' display_flag_text];
    savefile_loading_case = 'MultiAx';
else
    % X-Axis Bending
    eps0 = [0.1, 0, 0]';
    k0 = [0,0,0]';
    title_ = ['X-Axis Bending ($\epsilon_{0_{1}}=0.1$) - ' display_flag_text];
    savefile_loading_case = 'XBend';
end


filenames = ['ANALYSIS_DiagonalFlatColor_pk1';'ANALYSIS_DiagonalFlatColor_pk2'];


% Options to generate the combined Flat color picture
options.show_ticks = 0;
options.show_title = 1;
options.show_coords = {};
if cs_type == "circle"
    options.show_coords.flag = 0;
else
    options.show_coords.flag = 0;
end
options.show_coords.center = cs_coords_center;
options.given_title = "";%title_;
options.fontsize = 24;
options.A_text.pos = cs_A_text_pos;
options.A_text.text = "PK1";
options.A_text.font = 28;
options.A_text.color = "black"; %"white";
options.B_text.pos = cs_B_text_pos;
options.B_text.text = "PK2";
options.B_text.font = 28;
options.B_text.color = "black";%"white";
options.loading_case = loading_case;

if display_type == "vm"
    options.cb_decimals = 1;
    if cs_type == "circle"
        if loading_case == "full"
            options.cdata_mmin = 9.6;
            options.cdata_mmax = 14;
        else
            options.cdata_mmin = 0.76;
            options.cdata_mmax = 1.76;
        end
    else
        if loading_case == "full"
            options.cdata_mmin = 10.5;
            options.cdata_mmax = 13.6;
        else
            options.cdata_mmin = 1.1;
            options.cdata_mmax = 1.7;  
        end
    end
else
    options.cb_decimals = 3;
end

fnameA = [filenames(1, :), '.msh'];
fnameB = [filenames(2, :), '.msh'];





%%

dbc = [];        % dbc = [node index, node dof, prescribed displacement]
tbc = [];
tol = 1e-8;
dof = 3;
eltype = 30;

index_SVK_pk1 = 14; % Saint-Vernant Kirchhoff with PK1
index_SVK_pk2 = 114; % Saint-Vernant Kirchhoff with PK2

indexes = [index_SVK_pk1, index_SVK_pk2];
titles = ["SVK (PK1)", "SVK (PK2)", ...
    "Mooney-Rivlin (PK1)", "Mooney-Rivlin (PK2)", ...
    "Neo-Hook (PK1)", "Neo-Hook (PK2)"];
colors = ["blue", "blue", ...
    "green", "green", ...
    "red", "red"];
shapes = ["square", "square", "^", "^", "o", "o"];


if recompute_files == 1
    for j = 1:2
        mesh = build_iga_mesh( plate );
        curve = extract_iga_boundary(mesh);
        
        % Retrieve material properties
        mat = default_mat();
        mat.index = indexes(j);
        filename = filenames(j, :);
        fname = get_output_file_name(filename);
        fout = fopen(fname,'w'); 
        
        nl_returns = nliga_returns( eltype, plate, mesh, mat, dbc, tbc, fout, eps0, k0);
        disp("Done Bending " + num2str(j));
          
        fclose(fout);
    end
    pause(3); % Pauso to finish writing the output files to disk
end

plot_color_flat_combined(display_flag, fnameA, fnameB, options);

if save_file == 1
    savefile_name = join(['RESULT_' savefile_loading_case savefile_display_type savefile_cs]);
    savefile_path = fullfile(pwd, 'output', savefile_name)
    exportgraphics(gcf,savefile_path,'Resolution',300);
end
