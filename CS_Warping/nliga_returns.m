function [nliga_return] = nliga_returns( eltype, geo, mesh, mat, dbc, tbc, fout,eps0, k0 )
% Main frame for the nonlinear isogeometric analysis, expanded for
% compatiblity with the CSWP and supporting an output-format
% Input:
    % eltype    - Element type 
    %               10 - plane strain element
    %               20 - solid element
    %               30 - CSWP element
	% geo   - Employed IGA Geometry 
	% mesh  - Employed mesh 
	% mat   - (Struct) containing material parameters
    % dbc   - displacements boundary conditions
    % tbc   - tractions boundary conditions 
    % fout  - output figure handle for visualization
    % eps0  - (3,1) vector containing the strain prescriptors
    % k0    - (3,1) vector containing the twist prescriptors
% Output:
    % nliga_return  - (Struct) Containing different tracked variables and
%                       solution data
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

% initialization
if eltype == 10      % plane element
    dof = 2;    % degree of freedom
    egp = (mesh.p+1)*(mesh.q+1); % gauss point in each element
elseif eltype == 20  % solid element
    dof = 3;    % degree of freedom
    egp = (mesh.p+1)*(mesh.q+1)*(mesh.k+1); % gauss point in each element
elseif eltype == 30  % solid element
    dof = 3;    % degree of freedom
    egp = (mesh.p+1)*(mesh.q+1); % gauss point in each element
end
ndofs = dof * mesh.nCpts;   % total dofs
ngp   = egp*mesh.nElems;    % total integration points
scatdbc = [];
scattbc = [];
if ~isempty(dbc)
    scatdbc = dof * (dbc(:,1)-1) + dbc(:,2);   % scatter dbc
end
if ~isempty(tbc) 
    scattbc = dof * (tbc(:,1)-1) + tbc(:,2);   % scatter tbc
end
tol = 1e-6;           % tolerance of convergence
maxit = 20;           % maximum iterative steps 
reit = 0;             % reduction index
maxreit = 6;          % maximum load/displacement step reduction 
ndbc = size(dbc,1);   % number of displacement constrained nodes
ntbc = size(tbc,1);   % number of displacement constrained nodes
if eltype==30
    u = zeros(ndofs+6,1);
    cu = zeros(ndofs + 6,1);  % converged nodal displacements
else
    u = zeros(ndofs,1);   % nodal displacements
    cu = zeros(ndofs,1);  % converged nodal displacements
end


step = 0;             % load/displacement step index
curtime = 0;          % current time
timeInterval = 1;  % initial time interval
cnit = [];            % record the iterative steps

% belongs to hyperelastic materials
if ((mat.index >= 10 && mat.index < 20) || (mat.index >= 110 && mat.index < 120))            
    % % output initial undeformed geometries
    if eltype ==30 % CSWP-Element
        output_visual_mesh_CSWP( fout, mat, geo, mesh, u, step, curtime );
    elseif mesh.dim == 2 % Plane Element
        output_visual_mesh2d( fout, mat, geo, mesh, u, step, curtime );
    elseif mesh.dim == 3 % Block Element
        output_visual_mesh3d( fout, mat, geo, mesh, u, step, curtime );
    end
end

while curtime ~= 1    % get to the end
    curtime = curtime + timeInterval;
    if curtime > 1
        timeInterval = 1 - curtime + timeInterval;
        curtime = 1; 
    end
    err = 1e6;        % record iterative error
    perr = err;
    nit = 0;          % record iterative times
    fprintf(1,'\n \t time    time step    iter \t  residual \n');
    if eltype ==30
        iu = zeros(ndofs+6,1);   % record iterative displacements
    else
        iu = zeros(ndofs,1);   % record iterative displacements
    end
    while (err > tol) && ( nit <= maxit)
        nit = nit+1;
        % CSWP Element - belongs to the cross-sectional warping problem
        if eltype == 30 
            if mat.index >= 10 && mat.index < 20
                % Elastic CSWP with PK1
                [ k, r ] = globalstiffness_CSWP( eltype, geo, mesh, mat, u, curtime,eps0,k0 );
            elseif mat.index >= 110 && mat.index < 120
                % Elastic CSWP with PK2
                [ k, r ] = globalstiffness_CSWP_PK2( eltype, geo, mesh, mat, u , curtime,eps0,k0);
            end
        % belongs to hyperelastic materials
        elseif ( mat.index >= 10 && mat.index < 20 )            
            % ATTENTION: Depricated
            [ k, r ] = globalstiffness_hyper( eltype, geo, mesh, mat, u);
        % belongs to plastic materials
        elseif ( mat.index >= 20 && mat.index < 40 )        
            % ATTENTION: Depricated
            [ k, r ] = globalstiffness_plastic( D, eltype, geo, mesh, mat, iu );
        end
        if eltype == 30
            f = zeros(ndofs+6,1);         % define external force
        else
            f = zeros(ndofs,1);         % define external force
        end
        if ntbc~=0                  % enforce traction conditions  
            f(scattbc) = tbc(:,3);
        end
           
        if ndbc~=0                  % enforce displacement conditions
            k(scatdbc,:) = zeros(ndbc, ndofs);
            k(scatdbc,scatdbc) = eye(ndbc);
            f(scatdbc,:) = 0;
            if nit == 1
                f(scatdbc,:) = dbc(:,3); 
            end          
        end
        b = curtime*f - r;          % define right side of the governing equation
        if ndbc~=0  
            b(scatdbc) = curtime*dbc(:,3) - u(scatdbc);   
        end
        du = k\b;                   % solve equation
        
        alldof = 1:ndofs;
        freedof = setdiff(alldof, scatdbc);    % nodes without displacement constraint
        u = u + du;                 % update displacement 
        iu = iu + du;               % update increment displacement
        if nit > 1                  % compute iterative error
            num = b(freedof)' * b(freedof);
            denom = 1+f(freedof)' * f(freedof);
            err = num/denom; 
        end
        
        % output current time step and iterative error
        fprintf(1,'%10.5f %10.3e %5d %14.5e \n',curtime,timeInterval,nit,err); 
        if err/perr > 1E3 && nit > 2
            nit = maxit+1;   % lf solution diverge extremely, go to next iteration
        else
            perr = err;
        end
    end
    
    if  nit <= maxit               % converged 
        reit = 0;                  % reset reduction index
        step = step + 1;           % increase converged steps by 1
        cu = u;                    % update converged displacement
        cnit = [cnit, nit];
        if length(cnit) >=2 && all(cnit(end-1:end) <= 5)
            timeInterval = timeInterval*1.5;  % increase the increment by times 1.5 
        end
        % belongs to hyperelastic materials
        if ( mat.index >= 10 && mat.index < 20 ) || ( mat.index >= 110 && mat.index < 120 )            
            % output visualized mesh file with 'filename'
            %if mesh.dim == 2 && eltype == 30
            if eltype == 30
                output_visual_mesh_CSWP( fout, mat, geo, mesh, u, step, curtime );
            elseif mesh.dim == 2 
                output_visual_mesh2d( fout, mat, geo, mesh, u, step, curtime );
            elseif mesh.dim == 3
                output_visual_mesh3d( fout, mat, geo, mesh, u, step, curtime );
            end
        elseif ( mat.index >= 20 && mat.index < 40 )        % belongs to plastic materials   
            % Output for NLIGA-Plasticity Model
            output_plastic( D, eltype, fout, mat, geo, mesh, iu, u, step, curtime);
        end
        
    else                           % not converged
        if reit <= maxreit         % refine time interval and continue iterating
            curtime = curtime - timeInterval;   % recover current time step
            init_vina(ngp) ;                    % Reset the memory-variables ?
            timeInterval = timeInterval/4;      % refine time interval
            reit = reit+1;         % increase reduction index
            u = cu;                % recover current displacement from last converged displacement
        else
            return;                % stop analysis
        end
    end
end

% Determination of Beam Stiffness
C0 = beam_stiffness(geo, mesh, mat, eps0, k0, u, k);

    
% Determination of Beam Forces
[n0, m0] = beam_forces(geo, mesh, mat, eps0, k0, u);


% Determine the Deformed configuration
u_disp = u(1:end-6);
coords_def = mesh.coords(:, [1:3]) + reshape(u_disp, 3, [])';


% Store results in return data format
nliga_return.C0 = C0;
nliga_return.n0 = n0;
nliga_return.m0 = m0;
nliga_return.k = k;
nliga_return.u = u;
nliga_return.coords_def = coords_def;

end

