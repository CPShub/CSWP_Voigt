function mat = default_mat()
% Constructs a struct that holds basic material properties and
% the index marker (0 by default, changed as needed for different material models)
% Output:
	% mat   - (Struct) containing material index and properties
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

% Relevant material properties are:
%   lambda_mat:     121         % GPa
%   mue_mat:         90          % GPa
%   E_modulus:      208.1592    % GPa
%   Yield Stress:   250         % MPa
%   (Index):        0

mat.index = 0;
mat.lambda_mat = 121;
mat.mue_mat = 80; % GPa - Also called Shear Modulus

mat.E_modulus = mat.mue_mat * (3*mat.lambda_mat + 2*mat.mue_mat)/(mat.lambda_mat + mat.mue_mat);        % GPa
mat.poissons_ratio = (mat.lambda_mat)/ (2*mat.lambda_mat + 2*mat.mue_mat); % GPa
mat.yield_stress = 0.250;     % GPa

% Compression Modulus / Bulk Modulus
mat.compression_modulus = mat.E_modulus / (3*(1-2*mat.poissons_ratio));% 0 = Incompressible


% Neo-Hook Parameters
mat.NH_coef1 = mat.mue_mat / 2;

% Mooney-Rivlin Parameters
mat.MR_coef1 = 0.75 * mat.mue_mat / 2;
mat.MR_coef2 = 0.25 * mat.mue_mat / 2;

% YEOH-Parameters
mat.YEOH_coef1 = 1;
mat.YEOH_coef2 = 1;
mat.YEOH_coef3 = 1;

% BIDTANERMAN-Parameters
mat.BIDTANERMAN_coef1 = 1;
mat.BIDTANERMAN_coef2 = 1;
mat.BIDTANERMAN_coef3 = 1;
mat.BIDTANERMAN_coef4 = 1;
end