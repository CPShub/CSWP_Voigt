% This script is used to generate grey-scale images of the two used meshes: 
% unit square and unit circle
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


% Toggle to save the files directly into "output"
save_file = 1;





% O-Mesh Circle
cs_options = {};
cs_options.Refinement = 3;
cs_options.show_plot = 1;

plate = geo_circle_with_square( [0,0], 1, 0.6, cs_options);
axis off
xlim([-1.1, 1.1])
ylim([-1.1, 1.1])

% Assign the correct face color for the multi-patch
hSurfs = findobj(gcf, 'Type', 'surface');
delete(findobj(gcf, 'Type', 'light'));
set(hSurfs, 'FaceColor', 'flat');
set(hSurfs, 'FaceColor', "#aeaca4")

% Safe the File
if save_file == 1
    savefile_name = 'RESULT_CircleMesh.png';
    savefile_path = fullfile(pwd, 'output', savefile_name)
    exportgraphics(gcf,savefile_path,'Resolution',300);
end

% Square
cs_options = {};
cs_options.RefinementX = 9;
cs_options.RefinementY = 9;
cs_options.show_plot = 1;

figure
geo_square([0,0],1,cs_options)
axis off
xlim([-0.6, 0.6])
ylim([-0.6, 0.6])

% Safe the File
if save_file == 1
    savefile_name = 'RESULT_SquareMesh.png';
    savefile_path = fullfile(pwd, 'output', savefile_name)
    exportgraphics(gcf,savefile_path,'Resolution',300);
end
