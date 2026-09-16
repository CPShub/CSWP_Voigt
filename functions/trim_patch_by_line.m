function [trimA, trimB] = trim_patch_by_line(vmesh, ptA, ptB)
% Divides a vmesh structure along a line through points A and B
% Input:
	% vmesh     - output visualized mesh structure, see "read_visual_mesh()"
    % ptA       - (2,1) point in space
    % ptB       - (2,1) point in space
% Output:
    % trimA     - vmesh structure on the "positive" side of the diagonal
    % trimB     - vmesh structure on the "positive" side of the diagonal
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
% warping problem for hyperelastic beams: A compact formulation in 
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

d = ptB - ptA;
numSteps = length(vmesh.vertices);
e = 0;%1e-6;

trimA.vertices = cell(1, numSteps);
trimA.face = cell(1, numSteps);
trimA.displacement = cell(1, numSteps);
trimA.stress = cell(1, numSteps);
trimA.strain = cell(1, numSteps);

trimB = trimA;

for step = 1:numSteps
    % sub-select only the first num1*num2 points WITHOUT the kntcrv-points
    v_orig = vmesh.vertices{step};
    
    disp_orig = vmesh.displacement{step};
    stress_orig = vmesh.stress{step};
    strain_orig = vmesh.strain{step};

    v_rel = v_orig(:, 1:2) - ptA;
    sideValues = d(1) * v_rel(:, 2) - d(2) * v_rel(:, 1);

    if step == 1
        flagsA = true(size(sideValues));
        flagsB = true(size(sideValues));
    else
        flagsA = sideValues >= e;
        flagsB = sideValues <= -e;
    end
    
    trimA.vertices{step} = v_orig(flagsA, :);
    trimB.vertices{step} = v_orig(flagsB, :);

    trimA.displacement{step} = disp_orig(flagsA, :);
    trimB.displacement{step} = disp_orig(flagsB, :);
    
    trimA.stress{step} = stress_orig(flagsA, :);
    trimB.stress{step} = stress_orig(flagsB, :);
    
    trimA.strain{step} = strain_orig(flagsA, :);
    trimB.strain{step} = strain_orig(flagsB, :);
end
end