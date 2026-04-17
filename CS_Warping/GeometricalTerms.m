function [M,m,Theta] = GeometricalTerms(x0,x,mu,e)
% Calculate relevant geometrical terms for the cross-sectional warping
% Input:
	% x0    - Undeformed position vector
	% x     - Deformed position vector
	% mu    - (3,1) vector of lagrange multipliers (end-2:end of u)
    % e     - eye-matrix

% Output:
	% M     - Body Translation Tensor
	% m     - Body Rotation Tensor
    % Theta - Tensor Term
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

% The following expressions and equations can also be found in "Numerische 
% Methoden zur Modellierung elastoplastischer Balken und ihre Anwendung auf 
% periodische Gitterstrukturen", PhD Thesis by L. Herrnböck

% M = x(1)*(e(:,2)*e(:,3)'+e(:,3)*e(:,2)')+x(2)*(e(:,1)*e(:,3)'+e(:,3)*e(:,1)')+x(3)*(e(:,1)*e(:,2)'+e(:,2)*e(:,1)');
% m = x(2)*x(3)*e(:,1)+x(3)*x(1)*e(:,2)+x(1)*x(2)*e(:,3);
% Theta = mu(1)*(e(:,2)*e(:,3)'+e(:,3)*e(:,2)') + mu(2)*(e(:,1)*e(:,3)'+e(:,3)*e(:,1)') +mu(3)*(e(:,1)*e(:,2)'+e(:,2)*e(:,1)');

%Eq. 4.7
m = x(2)*x(3)*e(:,1)+x(3)*x(1)*e(:,2)+(atan(x(2)/x(1))-atan(x0(2)/x0(1)))*e(:,3);

%Eq. 4.10.
term0=x(1)^2+x(2)^2;
MT = [0 x(3) x(2); x(3) 0 x(1);-x(2)/term0 x(1)/term0 0];
M = MT';

%Eq. 4.16
term1=2*x(1)*x(2);
term2=x(2)^2-x(1)^2;
Theta = [mu(3)*term1/(term0^2) mu(3)*term2/(term0^2) mu(2); mu(3)*term2/(term0)^2 -mu(3)*term1/term0^2 mu(1); mu(2) mu(1) 0];
end