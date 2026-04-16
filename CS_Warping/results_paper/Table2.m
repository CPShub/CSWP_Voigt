% This script is used to compare the iterative errors in computation of the
% same CSWP-Problem using the PK1 and PK2 formulations. As is highlighted
% in the associated paper (1), both are (always) equal.
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

% Build geometrical model
plate =  geo_square( [0,0], 1, 0);
eps0 = [0.02, 0.03, 0.06]';
k0 = [0.01,0.02,0.1]';

% Store indicees for compared material models
index_SVK_pk1 = 14; % Saint-Vernant Kirchhoff with PK1
index_SVK_pk2 = 114; % Saint-Vernant Kirchhoff with PK2

indexes = [index_SVK_pk1, index_SVK_pk2];
l = length(indexes);

% Enforce displacement boundary conditions 
dbc = [];        
tbc = [];
dof = 3;
eltype = 30;

% Retrieve material properties
mat = default_mat();

for i = 1:l
    mat.index = indexes(i);
    
    filename = 'ANALYSIS_PK2_ErrorTable';
    fname = get_output_file_name(filename);
    fout = fopen(fname,'w'); 
        
    % Execute NLIGA with CSWP
    nl_returns = nliga_returns( eltype, plate, mesh, mat, dbc, tbc, fout, eps0, k0);

    % Comparing the Residual entries from both nliga-calls highlights the
    % equality in their results
end