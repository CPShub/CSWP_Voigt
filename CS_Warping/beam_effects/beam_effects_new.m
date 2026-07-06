function [n0, m0, C0] = beam_effects_new(geo, mesh, mat, eps0, k0, u, K)


% Integration over the whole domain
gp_x = mesh.p+1;        % number of integration points in x-direction
gp_y = mesh.q+1;        % number of integration points in y-direction
[gp, wgt] = gauss_quadrature(gp_x, gp_y);   % calculate integration points and their weights

dof = 3;
ndofs = dof * mesh.nCpts;      % total dofs

% Preallocate space
R_y_all = zeros(ndofs, 6);
stress_resultant = zeros(6,1);
e = eye(3);


dk0dq = [zeros(3,3), eye(3,3)];
deps0dq = [eye(3,3), zeros(3,3)];

% First Integration to Compute the Stress Resultants n0, m0 and the Displacement
% Sensitivites u,y
for el = 1:mesh.nElems                % loop over elements
    sctr = mesh.elNodeCnt(el,:);       % element control points index
    elDoma = mesh.elDoma(el,:);        % element parametric domain
    elCpts0 = mesh.initcoords(sctr,:); % initial coordinates of el cont points
    nn = length(sctr);                % number of control points for each element
    nnElem = nn*dof;                  % total dof for each element
    sctrB = zeros(1, nnElem);

    for i = 1:dof
        sctrB(i:dof:nnElem) = dof*(sctr-1) + i;  % displacement in i-th direction
    end
    
    elDisp = u(sctrB);
    elDisp = reshape(elDisp, dof, nn);
    elCpts(:,1:3)=elCpts0(:,1:3)+elDisp'; %here we actualize x+du
    
    for ipt = 1:size(gp,1)            % loop over integration points
        pt = gp(ipt,:);      % reference parametric coordinates for each integration point
        wt = wgt(ipt);       % weigths for each integration point
        gauPts = parameter_gauss_mapping( elDoma, pt );   % gauss integration mapping
        j1 = jacobian_gauss_mapping( elDoma );     % jacobian value for gauss mapping   
        [N,ders] = nurbs_derivatives( gauPts,geo, mesh );
        jmatrix = ders*elCpts0(:,1:dof-1); %Because the mapping is in 2D
        j2 = det(jmatrix);
        ders =  jmatrix \ ders;              
        fac = j1 *j2 * wt;  

        ders3D = zeros(3,size(elCpts,1));
        ders3D(1:2,:) = ders;

        x = N.*elCpts(:,1:dof)'; % Position of the current Gauss Point in the mesh
        x = sum(x,2);

        dx_alpha = elDisp * ders3D';
        F = def_gradient(eps0, k0, x, dx_alpha);

        if (det(F))<0 % Check determinate 
            warning('det(F) is negative, error')
        end
        
        % Compute the material response
        if (mat.index >= 10 && mat.index < 20) % Formulation with PK1
            error("This solution scheme for the beam effects and stiffness is only usable with a PK2 material formulation!")
        elseif (mat.index >= 110 && mat.index < 120) % Formulation with PK2
            [pk2, dtan] = material_CSWP_PK2_hyperelasticity(dof, mat, F);
        end

        % Compute BN
        BN = zeros(6, nn*3);
        
        % Indicees for blocks
        idx1 = 1:3:3*nn;
        idx2 = 2:3:3*nn;
        idx3 = 3:3:3*nn;

        C3 = [ (k0(3)*F(2,3)-k0(2)*F(3,3)), (k0(1)*F(3,3)-k0(3)*F(1,3)), (k0(2)*F(1,3)-k0(1)*F(2,3)) ];
        C5 = [ (k0(3)*F(2,2)-k0(2)*F(3,2)), (k0(1)*F(3,2)-k0(3)*F(1,2)), (k0(2)*F(1,2)-k0(1)*F(2,2)) ];
        C6 = [ (k0(3)*F(2,1)-k0(2)*F(3,1)), (k0(1)*F(3,1)-k0(3)*F(1,1)), (k0(2)*F(1,1)-k0(1)*F(2,1)) ];
        

        % Row 1
        BN(1, idx1) = F(1,1)*ders(1,:); 
        BN(1, idx2) = F(2,1)*ders(1,:); 
        BN(1, idx3) = F(3,1)*ders(1,:);
        
        % Row 2
        BN(2, idx1) = F(1,2)*ders(2,:); 
        BN(2, idx2) = F(2,2)*ders(2,:); 
        BN(2, idx3) = F(3,2)*ders(2,:);
        
        % Row 3
        BN(3, idx1) = C3(1)*N;          
        BN(3, idx2) = C3(2)*N;          
        BN(3, idx3) = C3(3)*N;
        
        % Row 4
        BN(4, idx1) = F(1,1)*ders(2,:) + F(1,2)*ders(1,:);
        BN(4, idx2) = F(2,1)*ders(2,:) + F(2,2)*ders(1,:);
        BN(4, idx3) = F(3,1)*ders(2,:) + F(3,2)*ders(1,:);
        
        % Row 5
        BN(5, idx1) = C5(1)*N + F(1,3)*ders(2,:);
        BN(5, idx2) = C5(2)*N + F(2,3)*ders(2,:);
        BN(5, idx3) = C5(3)*N + F(3,3)*ders(2,:);
        
        % Row 6
        BN(6, idx1) = F(1,3)*ders(1,:) + C6(1)*N;
        BN(6, idx2) = F(2,3)*ders(1,:) + C6(2)*N;
        BN(6, idx3) = F(3,3)*ders(1,:) + C6(3)*N;
        

        % Compute Derivatives of BN wrt q
        R_y_partial_all = zeros(nnElem, 6);
        
        % Extract the relevant stress components for PK2_reduced
        S_red = [pk2(6); pk2(5); pk2(3)];

        % Extract the relevant columns of the material stiffness matrix
        C_red = dtan(:, [3,5,6]);

        % compute the h,p operators for all variations of p in [eps0, k0]
        % Assemble them directly into the H operator
        H_eps0 = F';
        dkdq_cross_x = [0,x(3),-x(2); -x(3), 0, x(1); x(2), -x(1), 0];
        deps0k0_dq = [eye(3,3), dkdq_cross_x];
        H_k0 = F' * dkdq_cross_x;%cross(eye(3,3), repmat(x, 1,3)); %TODO: vectorize this!
        H = [H_eps0, H_k0];

        % Safe F column vectors
        F1 = F(:, 1);
        F2 = F(:, 2);
        F3 = F(:, 3);

        % Compute the E,q^{0} operator
        %dE3_dq = 
        dE_dq_zero = [zeros(3,6);
            H(3,:);
            zeros(1,6);
            H(2,:);
            H(1,:)];
        
        % % Compute the B_{i,q}^{0,T} operators for i = 3,5,6
        % BN_3q_0 = -N * ( [zeros(3,3), cross(eye(3,3), repmat(F3,1,3))] ...
        %     + cross(repmat(k0,1,6),deps0k0_dq));
        % 
        % BN_5q_0 = ders(2, :) * (deps0k0_dq) ...
        %     - N * [zeros(3,3), cross(eye(3,3), repmat(F2,1,3))];
        % 
        % BN_6q_0 = ders(1, :) * (deps0k0_dq) ...
        %     - N * [zeros(3,3), cross(eye(3,3), repmat(F1,1,3))];
        % 
        % BN_q_0T_red = zeros(3*nn, 3);

        %%%%%%%%%%%%%%%%%%%%%%%
        
        BN_q_0 = zeros(6, nn*3); % Ziel-Dimension (6 x 48)
        
        % Indizes für die Blöcke der 3 Freiheitsgrade (genau wie in deinem BN-Code)
        idx1 = 1:3:3*nn;
        idx2 = 2:3:3*nn;
        idx3 = 3:3:3*nn;
        
        for node = 1:nn
            % 1. Werte für den aktuellen Knoten isolieren
            N_i     = N(node);
            ders1_i = ders(1, node);
            ders2_i = ders(2, node);
            
            % 2. Die Kreuzprodukt-Ausdrücke für DIESEN Knoten berechnen
            % Da wir (6, 48) aufbauen, betrachten wir hier die ersten 3 Spalten (1:3) 
            % von deps0k0_dq, da diese zu den ersten 3 Freiheitsgraden gehören.
            
            % Für Zeile 3 (Dehnung 3):
            % [zeros(3,3), cross(...)] -> die ersten 3 Spalten sind 0.
            % Es bleibt für die ersten 3 Spalten nur der cross(k0, deps0k0_dq) Anteil.
            C3 = -N_i * ( cross(k0, deps0k0_dq(:, 1:3)) );
            
            % Für Zeile 5 (Dehnung 5):
            % Auch hier fällt der rechte zeros(3,3)-Anteil komplett weg!
            C5 = ders2_i * deps0k0_dq(:, 1:3);
            
            % Für Zeile 6 (Dehnung 6):
            % Auch hier fällt der rechte zeros(3,3)-Anteil weg.
            C6 = ders1_i * deps0k0_dq(:, 1:3);
            
            % 3. Zeile für Zeile in die globalen Spalten des Knotens einsortieren
            % Zeilen 1, 2 und 4 sind für diesen Anteil komplett 0 (da q_0 dort nicht eingeht)
            
            % Row 3
            BN_q_0(3, idx1(node)) = C3(1);          
            BN_q_0(3, idx2(node)) = C3(2);          
            BN_q_0(3, idx3(node)) = C3(3);
            
            % Row 5
            BN_q_0(5, idx1(node)) = C5(1);
            BN_q_0(5, idx2(node)) = C5(2);
            BN_q_0(5, idx3(node)) = C5(3);
            
            % Row 6
            BN_q_0(6, idx1(node)) = C6(1);
            BN_q_0(6, idx2(node)) = C6(2);
            BN_q_0(6, idx3(node)) = C6(3);
        end
        

        % for i = 1:nn
        %     % Corresponds to Eq. 52 in (1)
        %     BN(:,i*3-2:i*3) = [ 
        %     F(1,1)*ders(1,i)     F(2,1)*ders(1,i)      F(3,1)*ders(1,i);
        %     F(1,2)*ders(2,i)     F(2,2)*ders(2,i)      F(3,2)*ders(2,i);
        %     N(i)*(k0(3)*F(2,3)-k0(2)*F(3,3))    N(i)*(k0(1)*F(3,3)-k0(3)*F(1,3))       N(i)*(k0(2)*F(1,3)-k0(1)*F(2,3)) ;
        %     F(1,1)*ders(2,i)+ F(1,2)*ders(1,i)  F(2,1)*ders(2,i)+F(2,2)*ders(1,i)   F(3,1)*ders(2,i)+F(3,2)*ders(1,i);
        %     (F(1,3)*ders(2,i) + N(i)*(k0(3)*F(2,2)-k0(2)*F(3,2)))  (F(2,3)*ders(2,i) + N(i)*(k0(1)*F(3,2)-k0(3)*F(1,2)))   (F(3,3)*ders(2,i) + N(i)*(k0(2)*F(1,2)-k0(1)*F(2,2)));
        %     (F(1,3)*ders(1,i) + N(i)*(k0(3)*F(2,1)-k0(2)*F(3,1)))  (F(2,3)*ders(1,i) + N(i)*(k0(1)*F(3,1)-k0(3)*F(1,1)))   (F(3,3)*ders(1,i) + N(i)*(k0(2)*F(1,1)-k0(1)*F(2,1))) ];
        % end
        
        
        
        
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%




        % Compute R_y_partial_all 
        R_y_partial_all = BN' * C_red * H + BN_q_0T_red * S_red;

        % Combine with S_reduced and integrate
        stress_resultant = stress_resultant + fac * H' * S_red;

        % Integrate the sensitivities residual vector
        R_y_all(sctrB, :) = R_y_all(sctrB, :) + fac * R_y_partial_all;
    end

    % Extract n0 and m0
    n0 = stress_resultant(1:3)';
    m0 = stress_resultant(4:6)';
    C0 = zeros(6,6); %TODO: implemetnt this following eq. 132
end