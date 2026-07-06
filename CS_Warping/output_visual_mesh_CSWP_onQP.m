function vmesh = output_visual_mesh_CSWP_onQP( fout, mat, geo, mesh, u, step, currentime, eps0, k0)

% TODO: Write a DOC

% Check if mesh contains sub_meshes

if isfield(mesh, "submeshes")
    mesh_cell = mesh.submeshes;
    num_meshes = size(mesh.submeshes, 2);
    use_gloElNodeCnt = 1;
else
    mesh_cell = cell(1,1);
    mesh_cell{1,1} = mesh;
    num_meshes = 1;
    use_gloElNodeCnt = 0;
end

% Total number of points = Total number of gauss Quadrature Points
dof = 3;
gp_x = mesh.p+1;        % number of integration points in x-direction
gp_y = mesh.q+1;        % number of integration points in y-direction
[gp, ~] = gauss_quadrature(gp_x, gp_y);   % calculate integration points and their weights



% Pre-Assign vmesh
numpts = mesh.nElems * size(gp, 1);
vmesh.nodalpts = zeros(numpts,3);
vmesh.displacement = zeros(numpts,3);
vmesh.stress = zeros(numpts,6);
vmesh.strain = zeros(numpts,6);

% TODO: remove after testing
vms = zeros(numpts,1);

% Iterate over all meshes
j = 0;
for m = 1:num_meshes
    sub_mesh = mesh_cell{1, m};
    %sub_geo = geo{1,5};

    for e = 1:sub_mesh.nElems
        sctr = sub_mesh.elNodeCnt(e,:);     % element control points index
        %exyz = sub_mesh.coords(sctr,:);  % element control points' coordinates
        nn = length(sctr);   % number of control points in the element
        nn3 = nn*3;          % degree of freedom of control points
        %nn3 = nn*2;
        elDoma = sub_mesh.elDoma(e, :);
        
        
        % Check if global globElNodeCnt should be used
        if use_gloElNodeCnt
            sctr = sub_mesh.gloElNodeCnt(e, :);
        end
        sctrB = zeros(1, nn3);      
        sctrB(1:3:nn3) = 3*sctr - 2;% displacement in x direction
        sctrB(2:3:nn3) = 3*sctr-1;  % displacement in y direction
        sctrB(3:3:end) = 3*sctr;    % displacement in z direction
        edsp = u(sctrB);
        edsp = reshape(edsp, 3, nn);

        elCpts0 = mesh.initcoords(sctr,:); % initial coordinates of el cont points
        elCpts(:,1:3)=elCpts0(:,1:3)+edsp';

        for ipt = 1:size(gp,1)
            pt = gp(ipt,:);      % reference parametric coordinates for each integration point
            gauPts = parameter_gauss_mapping( elDoma, pt );   % gauss integration mapping  
            [N,ders] = nurbs_derivatives( gauPts,geo, mesh );
            jmatrix = ders*elCpts0(:,1:dof-1); %Because the mapping is in 2D
            ders =  jmatrix \ ders;    
            ders3D = zeros(3,size(elCpts,1));
            ders3D(1:2,:) = ders;
            x = N.*elCpts(:,1:dof)';
            x = sum(x,2);
    
            dx_alpha = edsp * ders3D';
            F = def_gradient(eps0, k0, x, dx_alpha);
    
            % Retrieve PK2 material response as PK2 Stress and dtangent 
            [ stress, ~ ] = material_CSWP_PK2_hyperelasticity( dof, mat, F );
            
            ccy = pk2cauchy(stress, F);
        
            % global vmesh Index j
            j = j + 1;
            vmesh.displacement(j,:) = (edsp * N')';
            vmesh.nodalpts(j,:) = x;%N*exyz(:,1:3) + vmesh.displacement(j,1:3);   
            vmesh.stress(j,:) = ccy';
            strain = (F'*F-eye(3))/2;
            vmesh.strain(j,:) = voigt(strain)';
    
            % TODO: remove after testing
            vms(j) = von_mises(ccy');
        end
    end
end

% Write to the output file
fprintf(fout,'STEP = %d, TIME = %e\n', step, currentime);
for i = 1:j
    fprintf(fout,'v %f %f %f\n', vmesh.nodalpts(i,:));
    fprintf(fout,'d %f %f %f\n', vmesh.displacement(i,:));
    fprintf(fout,'s %f %f %f %f %f %f\n', vmesh.stress(i,:));
    fprintf(fout,'t %f %f %f %f %f %f\n', vmesh.strain(i,:));
end
end

