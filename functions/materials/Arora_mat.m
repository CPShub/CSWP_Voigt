function mat = Arora_mat()
% Constructs a struct that holds basic material properties and
% the index marker (0 by default, changed as needed for different material models)
% This Material is adapted to be identical to the one used by Arora (see
% (1) for more information)
% Output:
	% mat   - (Struct) containing material index and bais properties
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

% Relevant material properties are:
%   lambda_mat:     121         % GPa
%   mu_mat:         80.194      % GPa
%   E_modulus:      206.768     % GPa
%   shear_modulus:  0.3010      % Gpa
%   Yield Stress:   450         % Gpa
%   (Index):        0

mat.index = 0;

mat.mue_mat = 1.0;
mat.lambda_mat = 1.275;
mat.E_modulus = mat.mue_mat * (3*mat.lambda_mat + 2*mat.mue_mat) / (mat.lambda_mat + mat.mue_mat);
mat.poissons_ratio = mat.lambda_mat / (2 * (mat.lambda_mat + mat.mue_mat));
mat.compression_modulus = mat.lambda_mat + (2/3) * mat.mue_mat;

% Plastic orthropic factors
mat.plastic_mats = [0.672790077, 0.672790077, 0.666666667, 0, 0, 0, 0.785052599, 0.785052599, 1.012246821];

% Material hardening parameters (Example values)
mat.k = 0.01;
mat.h = 0.002;

% Plastic spin parameter
mat.eta_spin = -500;

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