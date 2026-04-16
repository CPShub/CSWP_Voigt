function vmesh = output_visual_mesh_CSWP( fout, mat, geo, mesh, u, step, currentime )
% This function generates the virtual mesh (vmesh) struct and writes its
% contents to a .msh file. 
% This variant is meant to work with results generated solving the CSWP
% Input:
    % fout          - file handle of visualized mesh
	% mat           - (Struct) containing material parameters
	% geo           - Employed IGA Geometry 
	% mesh          - Employed mesh 
    % u             - Displacement solution vector
    % step          - Current simulation step
    % currentime    - Current simulation timestamp   
% Output:
	% vmesh         - output visualized mesh structure, see "read_visual_mesh()"
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

num1 = 40;%100;        
num2 = 40;%10;
polygon = build_visual_mesh_suf( num1, num2 );        % build visualized mesh
kntcrv = build_visual_knotcurve_suf( mesh.uKnots, mesh.vKnots, num1+1 ); % build visualized knot curves

numpts = (num1+1)*(num2+1) + size(kntcrv.linpts,1);
vmesh.nodalpts = zeros(numpts,2);
vmesh.displacement = zeros(numpts,3);
vmesh.stress = zeros(numpts,3);
vmesh.strain = zeros(numpts,3);
elem_index = find_point_span( mesh, polygon.tripts );
for i=1:(num1+1)*(num2+1)
    xi = polygon.tripts(i,1);  % u mesh point
    eta = polygon.tripts(i,2);  % v mesh point
    e = elem_index(i);   % element number
    sctr = mesh.elNodeCnt(e,:);     % element control points index
    exyz = mesh.coords(sctr,:);  % element control points' coordinates
    nn = length(sctr);   % number of control points in the element
    nn3 = nn*3;          % degree of freedom of control points
    %nn3 = nn*2;
    sctrB = zeros(1, nn3);      
    sctrB(1:3:nn3) = 3*sctr - 2;% displacement in x direction
    sctrB(2:3:nn3) = 3*sctr-1;  % displacement in y direction
    sctrB(3:3:end) = 3*sctr;    % displacement in z direction
    edsp = u(sctrB);
    edsp = reshape(edsp, 3, nn);

    [R,dRdparm] = nurbs_derivatives( [xi, eta],geo, mesh );
    jmatrix = dRdparm*exyz(:,1:2); 
    dRdx =  jmatrix \ dRdparm;              
    f = edsp(1:2,:) * dRdx' + eye(2);   

    [ pk2, ~ ] = constitutive_relation( mesh.dim, mat, f );
    ccy = pk2cauchy(pk2, f);
    vmesh.displacement(i,:) = edsp * R';
    vmesh.nodalpts(i,:) = R*exyz(:,1:2) + vmesh.displacement(i,1:2);   
    vmesh.stress(i,:) = ccy;
    strain = (f'*f-eye(2))/2;
    vmesh.strain(i,:) = voigt(strain);
end

count = (num1+1)*(num2+1);
line_index = find_point_span( mesh, kntcrv.linpts );
kntcrv.linmesh = kntcrv.linmesh + (num1+1)*(num2+1);
for i = 1:size(kntcrv.linpts,1)
    count = count+1;
    xi = kntcrv.linpts(i,1);
    eta = kntcrv.linpts(i,2);
    e = line_index(i);   % element number
    sctr = mesh.elNodeCnt(e,:);     % element control points index
    exyz = mesh.coords(sctr,:);  % element control points' coordinates
    nn = length(sctr);   % number of control points in the element
    nn3 = nn*3;          % degree of freedom of control points %% Check the 3 out!!! replace for 2
    %nn3 = nn*2;
    sctrB = zeros(1, nn3); 
    sctrB(1:3:nn3) = 3*sctr - 2;% displacement in x direction
    sctrB(2:3:nn3) = 3*sctr-1;  % displacement in y direction
    sctrB(3:3:end) = 3*sctr;    % displacement in z direction
    edsp = u(sctrB);
    edsp = reshape(edsp, 3, nn);

    [R,dRdparm] = nurbs_derivatives( [xi, eta],geo, mesh );
    jmatrix = dRdparm*exyz(:,1:2); 
    dRdx =  jmatrix \ dRdparm;              
    f = edsp(1:2, :) * dRdx' + eye(2);   
    [ pk2, ~ ] = constitutive_relation( mesh.dim, mat, f );
    cauchy = pk2cauchy( pk2, f );
    vmesh.displacement(count,:) = edsp * R';
    vmesh.nodalpts(count,:) = R*exyz(:,1:2) + vmesh.displacement(count,1:2);   
    vmesh.stress(count,:) = cauchy;
    strain = (f'*f-eye(2))/2;
    vmesh.strain(i,:) = voigt(strain);
end

fprintf(fout,'STEP = %d, TIME = %e\n', step, currentime);
for i = 1:numpts
    fprintf(fout,'v %f %f\n', vmesh.nodalpts(i,:));
    fprintf(fout,'d %f %f %f\n', vmesh.displacement(i,:));
    fprintf(fout,'s %f %f %f\n', vmesh.stress(i,:));
    fprintf(fout,'t %f %f %f\n', vmesh.strain(i,:));
end
for i = 1:size(polygon.trimesh,1)
    fprintf(fout,'f %d %d %d\n', polygon.trimesh(i,:));
end
for i = 1:size(kntcrv.linmesh,1)
    fprintf(fout, 'l ');
    for j = 1:size(kntcrv.linmesh,2)
        fprintf(fout,'%d ', kntcrv.linmesh(i,j));
    end
    fprintf(fout, '\n');
end


end

