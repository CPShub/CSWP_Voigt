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

% 1. Compute the Beam Forces acting on the cros-section
[forces, moments] = beam_forces(geo, mesh, mat, eps0, k0, u);
disp("Forces in [x,y,z]: ")
disp(forces);
disp("Moments in [x,y,z]: ")
disp(moments);


% 2. Compute the Beam Stiffness Matrix (Sensitivity of Forces and Moments in relation
% to eps0 and k0)
[C0] = beam_stiffness(geo, mesh, mat, eps0, k0, u, k);
disp("Beam Stiffness Matrix [6,6]: ")
disp(C0);
