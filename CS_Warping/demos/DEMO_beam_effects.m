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


% Solve the NLIGA simulation
nl_return = nliga_returns(eltype, geo, mesh, mat, dbc, tbc, fout, eps0, k0);
u = nl_return.u; % Solution displacement vector
k = nl_return.k; % Solution stiffness matrix

% Compute the Beam Forces and Moments acting on the cross-section as well 
% as the deformation solution sensitivities u,q and the Beam Stiffness Matrix 
[forces, moments, stiffness, sensitivities] = beam_effects(geo, mesh, mat, eps0, k0, u, k);

disp("Forces in [x,y,z]: ")
disp(forces);
disp("Moments in [x,y,z]: ")
disp(moments);
disp("Beam Stiffness Matrix [6,6]: ")
disp(stiffness);
