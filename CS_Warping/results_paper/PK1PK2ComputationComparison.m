% This script is used to compare the computational results and demands
% between the Arora-based PK1 formulation and the novel PK2 formulation.
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


%% Options

plot_K_matricees = 1;
plot_mesh_computation_time = 0;
plot_K_computation_comparison = 0;
plot_K_computation_comparison_single = 0;
%%



%% Visualization

if plot_K_matricees

    % Build geometrical model
    plate =  geo_square( [0,0], 1, 0);
    mesh = build_iga_mesh( plate );
    eps0 = [0.02, 0.03, 0.06]';
    k0 = [0.01,0.02,0.1]';
    
    % Store indicees for compared material models
    index_SVK_pk1 = 14; % Saint-Venant Kirchhoff with PK1
    index_SVK_pk2 = 114; % Saint-Venant Kirchhoff with PK2
    
    % Enforce displacement boundary conditions 
    dbc = [];        
    tbc = [];
    dof = 3;
    eltype = 30;
    
    % Retrieve material properties
    mat = default_mat();
    mat_pk1 = mat;
    mat_pk1.index = index_SVK_pk1;
    mat_pk2 = mat;
    mat_pk2.index = index_SVK_pk2;
    
    filename_pk1 = 'ANALYSIS_PK1_ComputationComparison';
    filename_pk2 = 'ANALYSIS_PK2_ComputationComparison';
    
    fout_pk1 = fopen(get_output_file_name(filename_pk1),'w');
    fout_pk2 = fopen(get_output_file_name(filename_pk2), 'w');
        
        
    % Execute NLIGA with CSWP
    nl_returns_pk1 = nliga_returns( eltype, plate, mesh, mat_pk1, dbc, tbc, fout_pk1, eps0, k0);
    nl_returns_pk2 = nliga_returns( eltype, plate, mesh, mat_pk2, dbc, tbc, fout_pk2, eps0, k0);
    
    K_pk1 = nl_returns_pk1.k;
    K_pk2 = nl_returns_pk2.k;



    fig = figure("Color","w", "Position", [100, 100, 1200, 400]);

    ticks_vals = 0:20:120;
    font_name = 'Helvetica';
    font_size = 10;

    % Subplot 1
    subplot(1,3,1)
    imagesc(K_pk1)
    axis square
    ax = gca;
    ax.YDir = 'reverse';
    ax.XAxisLocation = 'top';
    ax.XTick = ticks_vals;
    ax.YTick = ticks_vals;
    ax.FontName = font_name;
    ax.FontSize = font_size;
    ax.LineWidth = 1.0;
    title("Stiffness Matrix K_{PK1}", 'FontName', font_name, 'FontSize', font_size+2, 'FontWeight', 'bold')
    colormap(gca, bluewhitered(256))
    max_val1 = max(abs(K_pk1(:)));
    %if max_val1 > 0, clim([-max_val1, max_val1]); end
    cb1 = colorbar();
    cb1.FontName = font_name;
    cb1.FontSize = font_size;

    % Subplot 2
    subplot(1,3,2)
    diff_K = K_pk2 - K_pk1;
    imagesc(diff_K)
    axis square
    ax = gca;
    ax.YDir = 'reverse';
    ax.XAxisLocation = 'top';
    ax.XTick = ticks_vals;
    ax.YTick = ticks_vals;
    ax.FontName = font_name;
    ax.FontSize = font_size;
    ax.LineWidth = 1.0;
    title("Difference K_{PK2} - K_{PK1}", 'FontName', font_name, 'FontSize', font_size+2, 'FontWeight', 'bold')
    colormap(gca, bluewhitered(256))
    max_val_diff = max(abs(diff_K(:)));
    %if max_val_diff > 0, clim([-max_val_diff, max_val_diff]); end
    cb2 = colorbar();
    cb2.FontName = font_name;
    cb2.FontSize = font_size;

    % Subplot 3
    subplot(1,3,3)
    imagesc(K_pk2)
    axis square
    ax = gca;
    ax.YDir = 'reverse';
    ax.XAxisLocation = 'top';
    ax.XTick = ticks_vals;
    ax.YTick = ticks_vals;
    ax.FontName = font_name;
    ax.FontSize = font_size;
    ax.LineWidth = 1.0;
    title("Stiffness Matrix K_{PK2}", 'FontName', font_name, 'FontSize', font_size+2, 'FontWeight', 'bold')
    colormap(gca, bluewhitered(256))
    max_val2 = max(abs(K_pk2(:)));
    %if max_val2 > 0, clim([-max_val2, max_val2]); end
    cb3 = colorbar();
    cb3.FontName = font_name;
    cb3.FontSize = font_size;
end








if plot_mesh_computation_time
    max_refinement = 5; 
    min_refinement = 0;
        
    refinements = min_refinement:max_refinement;
    l = length(refinements);
    max_iterations = 10; 
    n_elements = zeros(1,l);
    
    time_pk1_max_val = zeros(1,l);
    time_pk1_min_val = zeros(1,l);
    time_pk1_mean_val = zeros(1,l);
    time_pk1_std_val = zeros(1,l);
    time_pk2_max_val = zeros(1,l);
    time_pk2_min_val = zeros(1,l);
    time_pk2_mean_val = zeros(1,l);
    time_pk2_std_val = zeros(1,l);
    
    for i = 1:l
        ref = refinements(i);
        geo_options.RefinementX = ref;
        geo_options.RefinementY = ref;
        plate = geo_square([0,0], 1, geo_options);
        mesh = build_iga_mesh(plate);
        n_elements(i) = mesh.nElems;
        times_pk1 = zeros(1, max_iterations);
        times_pk2 = zeros(1, max_iterations);
        for j = 1:max_iterations
            tic_pk1 = tic;
            nl_returns_pk1 = nliga_returns(eltype, plate, mesh, mat_pk1, dbc, tbc, fout_pk1, eps0, k0);
            time_pk1 = toc(tic_pk1);
            times_pk1(j) = time_pk1;
            tic_pk2 = tic;
            nl_returns_pk2 = nliga_returns(eltype, plate, mesh, mat_pk2, dbc, tbc, fout_pk2, eps0, k0);
            time_pk2 = toc(tic_pk2);
            times_pk2(j) = time_pk2;
        end
        time_pk1_max_val(i) = max(times_pk1);
        time_pk1_min_val(i) = min(times_pk1);
        time_pk1_mean_val(i) = mean(times_pk1);
        time_pk1_std_val(i) = std(times_pk1);
        time_pk2_max_val(i) = max(times_pk2);
        time_pk2_min_val(i) = min(times_pk2);
        time_pk2_mean_val(i) = mean(times_pk2);
        time_pk2_std_val(i) = std(times_pk2);
    end

    fig = figure("Color","w", "Position", [100, 100, 700, 500]);
    font_name = 'Helvetica';
    font_size = 11;
    
    hold on
    grid on
    set(gca, 'XScale', 'log', 'YScale', 'log', 'FontName', font_name, 'FontSize', font_size, 'LineWidth', 1.0);
    
    x_fill = [n_elements, fliplr(n_elements)];
    
    % PK1 Ranges
    pk1_min_max_y = [time_pk1_min_val, fliplr(time_pk1_max_val)];
    pk1_std_y = [(time_pk1_mean_val - time_pk1_std_val), fliplr(time_pk1_mean_val + time_pk1_std_val)];
    
    fill(x_fill, pk1_min_max_y, [0.9, 0.9, 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    fill(x_fill, pk1_std_y, [0.7, 0.7, 0.7], 'EdgeColor', 'none', 'FaceAlpha', 0.6);
    
    % PK2 Ranges
    pk2_min_max_y = [time_pk2_min_val, fliplr(time_pk2_max_val)];
    pk2_std_y = [(time_pk2_mean_val - time_pk2_std_val), fliplr(time_pk2_mean_val + time_pk2_std_val)];
    
    fill(x_fill, pk2_min_max_y, [0.85, 0.9, 0.95], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    fill(x_fill, pk2_std_y, [0.65, 0.75, 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    
    % Means and Markers
    p1 = plot(n_elements, time_pk1_mean_val, '--k', 'LineWidth', 1.8);
    scatter(n_elements, time_pk1_mean_val, 40, 'k', 'filled');
    
    p2 = plot(n_elements, time_pk2_mean_val, '-b', 'LineWidth', 1.8);
    scatter(n_elements, time_pk2_mean_val, 40, 'b', 'filled');
    
    xlabel('Number of Elements, N', 'FontName', font_name, 'FontSize', font_size+1, 'FontWeight', 'bold');
    ylabel('Computation Time [s]', 'FontName', font_name, 'FontSize', font_size+1, 'FontWeight', 'bold');
    title('Computational Demand Comparison', 'FontName', font_name, 'FontSize', font_size+3, 'FontWeight', 'bold');
    
    legend([p1, p2], {'Arora-based PK1 formulation (Mean)', 'Novel PK2 formulation (Mean)'}, ...
        'Location', 'northwest', 'FontName', font_name, 'FontSize', font_size);
        
    hold off
end









%% plot_K_computation_comparison_single

if plot_K_computation_comparison_single
    % Setup all required information for the two computations
    max_iterations = 1000;

    % Material Properties
    index_SVK_pk1 = 14; % Saint-Venant Kirchhoff with PK1
    index_SVK_pk2 = 114; % Saint-Venant Kirchhoff with PK2
    mat = default_mat();
    mat_pk1 = mat;
    mat_pk1.index = index_SVK_pk1;
    mat_pk2 = mat;
    mat_pk2.index = index_SVK_pk2;
    
    geo_options.RefinementX = 2;
    geo_options.RefinementY = 2;
    geo_options.Degelev = [2, 2];
    plate = geo_square([0,0], 1, geo_options);
    mesh = build_iga_mesh(plate);
    
    gp_x = mesh.p+1;        % number of integration points in x-direction
    gp_y = mesh.q+1;        % number of integration points in y-direction
    [gp, wgt] = gauss_quadrature(gp_x, gp_y);   % calculate integration points and their weights
    
    %n_elements(i) = mesh.nElems;
    nn_gp = 1;%size(gp, 1);
    nn_elems = 1;%mesh.nElems;
    l = nn_elems * nn_gp;
    i = 0;

    % Time-Parameters for the Core
    times_pk1 = {};
    times_pk2 = {};
    times_pk1(l).mat = {};
    times_pk2(l).core ={};

    eps0 = [0.02, 0.03, 0.06]';
    k0 = [0.01,0.02,0.1]';
    
    dof=3;
    k0x = [    0    -k0(3)   k0(2);
        k0(3)      0    -k0(1);
       -k0(2)   k0(1)      0  ];
    e = eye(3,3);

    for el = 1:nn_elems
        sctr = mesh.elNodeCnt(el,:);
        elDoma = mesh.elDoma(el,:);        % element parametric domain
        elCpts0 = mesh.initcoords(sctr,:); % initial coordinates of el cont points
        nn = length(sctr);                % number of control points for each element
        nnElem = nn*dof;                  % dof for each element
        
        lambda = zeros(3,1);
        mu = zeros(3,1);
        elCpts(:,1:3)=elCpts0(:,1:3); %here we actualize x+du
        for ipt = 1:nn_gp
            i = i + 1;
            pt = gp(ipt,:);      % reference parametric coordinates for each integration point
            gauPts = parameter_gauss_mapping( elDoma, pt );   % gauss integration mapping  
            [N,ders] = nurbs_derivatives( gauPts,plate, mesh );
            jmatrix = ders*elCpts0(:,1:dof-1); %Because the mapping is in 2D
            ders =  jmatrix \ ders;                     
            ders3D = zeros(3,size(elCpts,1));
            ders3D(1:2,:) = ders;
            x = N.*elCpts(:,1:dof)';
            x = sum(x,2);
            x0 = N.*elCpts0(:,1:dof)';
            x0 = sum(x0,2);
            [M,m,Theta] = GeometricalTerms(x0,x,mu,e);
        
            dx_alpha = zeros(3,3);%elDisp * ders3D';
            F = def_gradient(eps0, k0, x, dx_alpha);
           
            times_pk2_core = zeros(max_iterations, 1);
            times_pk1_core = zeros(max_iterations, 1);
            times_pk2_mat = zeros(max_iterations, 1);
            times_pk1_mat = zeros(max_iterations, 1);
            for j = 1:max_iterations
                % Compute PK1 Material response
                tic
                [pk1, A] = material_CSWP_hyperelasticity(dof, mat_pk1, F);
                times_pk1_mat(j) = toc(tic);
            
                % Compute PK2 Material response
                tic
                [stress, dtan ] = material_CSWP_PK2_hyperelasticity( dof, mat_pk2, F );
                times_pk2_mat(j) = toc(tic);
        
                % Compute Core Responses
                [EK_pk2, R_pk2, time_pk2] = core_pk2(k0, F, stress, dtan, nnElem, nn, N, ders, k0x, lambda, mu, M, e, Theta);
                [EK_pk1, R_pk1, time_pk1] = core_pk1(pk1, A, nnElem, nn, N, ders, k0x, lambda, mu, M, e, Theta);
                times_pk2_core(j) = time_pk2;
                times_pk1_core(j) = time_pk1;
            end

            times_pk1(i).mat = times_pk1_mat;
            times_pk1(i).core = times_pk1_core;

            times_pk2(i).mat = times_pk2_mat;
            times_pk2(i).core = times_pk2_core;
        end
    end
    % Extract the arrays for the single-case (i = 1)
    pk1_core_data = times_pk1_core;%times_pk1(1).core;
    pk2_core_data = times_pk2_core;%times_pk2(1).core;
    pk1_mat_data  = times_pk1_mat;%times_pk1(1).mat;
    pk2_mat_data  = times_pk2_mat;%times_pk2(1).mat;

    fig = figure("Color", "w", "Position", [100, 100, 700, 700]);
    font_name = 'Helvetica';
    font_size = 11;
    
    pk1_fill_color = [0.90, 0.40, 0.40]; 
    pk1_line_color = [0.75, 0.15, 0.15]; 
    pk2_fill_color = [0.40, 0.60, 0.90]; 
    pk2_line_color = [0.15, 0.35, 0.75]; 
    
    cutoff = 1.1e-6;
    num_bins = 11; 
    shared_bins = linspace(0, cutoff, num_bins + 1);
    
    subplot(2,1,1)
    hold on
    grid on
    set(gca, 'FontName', font_name, 'FontSize', font_size, 'LineWidth', 1.0);
    
    pk1_core_clamped = min(pk1_core_data, cutoff);
    pk2_core_clamped = min(pk2_core_data, cutoff);
    
    h1 = histogram(pk1_core_clamped, 'BinEdges', shared_bins, 'Normalization', 'percentage', ...
        'FaceColor', pk1_fill_color, 'EdgeColor', pk1_line_color, 'FaceAlpha', 0.4, 'LineWidth', 1.0);
    
    h2 = histogram(pk2_core_clamped, 'BinEdges', shared_bins, 'Normalization', 'percentage',  ...
        'FaceColor', pk2_fill_color, 'EdgeColor', pk2_line_color, 'FaceAlpha', 0.4, 'LineWidth', 1.0);
    
    xlim([0, cutoff]);
    ylims_core = ylim;
    
    mean_pk1_core = mean(pk1_core_data);
    mean_pk2_core = mean(pk2_core_data);
    
    m1 = plot([mean_pk1_core, mean_pk1_core], ylims_core, '--', 'Color', pk1_line_color, 'LineWidth', 2.5);
    m2 = plot([mean_pk2_core, mean_pk2_core], ylims_core, '--', 'Color', pk2_line_color, 'LineWidth', 2.5);
    
    xlabel('Core Computation Time [s]', 'FontName', font_name, 'FontSize', font_size+1, 'FontWeight', 'bold');
    ylabel('Occurance in [%]', 'FontName', font_name, 'FontSize', font_size+1, 'FontWeight', 'bold');
    title('Core Computation Time Distribution', 'FontName', font_name, 'FontSize', font_size+2, 'FontWeight', 'bold');
    legend([h1, h2, m1, m2], {'PK1 Formulation (Core)', 'Novel PK2 Formulation (Core)', 'Mean PK1', 'Mean PK2'}, 'Location', 'northwest');
    hold off
    
    subplot(2,1,2)
    hold on
    grid on
    set(gca, 'FontName', font_name, 'FontSize', font_size, 'LineWidth', 1.0);
    
    pk1_mat_clamped = min(pk1_mat_data, cutoff);
    pk2_mat_clamped = min(pk2_mat_data, cutoff);
    
    h3 = histogram(pk1_mat_clamped, 'BinEdges', shared_bins, 'Normalization', 'percentage', ...
        'FaceColor', pk1_fill_color, 'EdgeColor', pk1_line_color, 'FaceAlpha', 0.4, 'LineWidth', 1.0);
    
    h4 = histogram(pk2_mat_clamped, 'BinEdges', shared_bins, 'Normalization', 'percentage', ...
        'FaceColor', pk2_fill_color, 'EdgeColor', pk2_line_color, 'FaceAlpha', 0.4, 'LineWidth', 1.0);
    
    xlim([0, cutoff]);
    ylims_mat = ylim;
    
    mean_pk1_mat = mean(pk1_mat_data);
    mean_pk2_mat = mean(pk2_mat_data);
    
    m3 = plot([mean_pk1_mat, mean_pk1_mat], ylims_mat, '--', 'Color', pk1_line_color, 'LineWidth', 2.5);
    m4 = plot([mean_pk2_mat, mean_pk2_mat], ylims_mat, '--', 'Color', pk2_line_color, 'LineWidth', 2.5);
    
    xlabel('Material Response Computation Time [s]', 'FontName', font_name, 'FontSize', font_size+1, 'FontWeight', 'bold');
    ylabel('Occurance in [%]', 'FontName', font_name, 'FontSize', font_size+1, 'FontWeight', 'bold');
    title('Material Response Computation Time Distribution', 'FontName', font_name, 'FontSize', font_size+2, 'FontWeight', 'bold');
    legend([h3, h4, m3, m4], {'PK1 Formulation (Mat)', 'Novel PK2 Formulation (Mat)', 'Mean PK1', 'Mean PK2'}, 'Location', 'northwest');
    hold off
end



if plot_K_computation_comparison_distro
    % Setup all required information for the two computations
    
    max_ipt = 5;
    min_ipt = 0;
        
    %refinements = min_refinement:max_refinement;
    %degelevs = min_degelev:max_degelev;
    ipt_elevs = min_ipt:max_ipt;
    l = length(ipt_elevs);
    max_iterations = 100; 
    %n_elements = zeros(1,l);

    % Material Properties
    index_SVK_pk1 = 14; % Saint-Venant Kirchhoff with PK1
    index_SVK_pk2 = 114; % Saint-Venant Kirchhoff with PK2
    mat = default_mat();
    mat_pk1 = mat;
    mat_pk1.index = index_SVK_pk1;
    mat_pk2 = mat;
    mat_pk2.index = index_SVK_pk2;
    
    
    % Time-Parameters for the Core
    time_pk1_max_val = zeros(1,l);
    time_pk1_min_val = zeros(1,l);
    time_pk1_mean_val = zeros(1,l);
    time_pk1_std_val = zeros(1,l);
    
    time_pk2_max_val = zeros(1,l);
    time_pk2_min_val = zeros(1,l);
    time_pk2_mean_val = zeros(1,l);
    time_pk2_std_val = zeros(1,l);

    % Time-Parameters for the Mat Response
    time_pk1_mat_max_val = zeros(1,l);
    time_pk1_mat_min_val = zeros(1,l);
    time_pk1_mat_mean_val = zeros(1,l);
    time_pk1_mat_std_val = zeros(1,l);
    
    time_pk2_mat_max_val = zeros(1,l);
    time_pk2_mat_min_val = zeros(1,l);
    time_pk2_mat_mean_val = zeros(1,l);
    time_pk2_mat_std_val = zeros(1,l);
    
    for i = 1:l
        val = ipt_elevs(i);
        geo_options.RefinementX = 2;
        geo_options.RefinementY = 2;
        geo_options.Degelev = [2, 2];
        plate = geo_square([0,0], 1, geo_options);
        mesh = build_iga_mesh(plate);
        
        gp_x = mesh.p+val;        % number of integration points in x-direction
        gp_y = mesh.q+val;        % number of integration points in y-direction
        [gp, wgt] = gauss_quadrature(gp_x, gp_y);   % calculate integration points and their weights
        
        %n_elements(i) = mesh.nElems;
        
        eps0 = [0.02, 0.03, 0.06]';
        k0 = [0.01,0.02,0.1]';
        
        dof=3;
        k0x = [    0    -k0(3)   k0(2);
            k0(3)      0    -k0(1);
           -k0(2)   k0(1)      0  ];
        e = eye(3,3);
        sctr = mesh.elNodeCnt(1,:);
        elDoma = mesh.elDoma(1,:);        % element parametric domain
        elCpts0 = mesh.initcoords(sctr,:); % initial coordinates of el cont points
        nn = length(sctr);                % number of control points for each element
        nnElem = nn*dof;                  % dof for each element
        
        lambda = zeros(3,1);
        mu = zeros(3,1);
        elCpts(:,1:3)=elCpts0(:,1:3); %here we actualize x+du
        ipt=1;
        pt = gp(ipt,:);      % reference parametric coordinates for each integration point
        gauPts = parameter_gauss_mapping( elDoma, pt );   % gauss integration mapping  
        [N,ders] = nurbs_derivatives( gauPts,plate, mesh );
        jmatrix = ders*elCpts0(:,1:dof-1); %Because the mapping is in 2D
        ders =  jmatrix \ ders;                     
        ders3D = zeros(3,size(elCpts,1));
        ders3D(1:2,:) = ders;
        x = N.*elCpts(:,1:dof)';
        x = sum(x,2);
        x0 = N.*elCpts0(:,1:dof)';
        x0 = sum(x0,2);
        [M,m,Theta] = GeometricalTerms(x0,x,mu,e);
    
        dx_alpha = zeros(3,3);%elDisp * ders3D';
        F = def_gradient(eps0, k0, x, dx_alpha);
       
        times_pk2_core = zeros(max_iterations, 1);
        times_pk1_core = zeros(max_iterations, 1);
        times_pk2_mat = zeros(max_iterations, 1);
        times_pk1_mat = zeros(max_iterations, 1);
        for j = 1:max_iterations
            % Compute PK1 Material response
            tic
            [pk1, A] = material_CSWP_hyperelasticity(dof, mat_pk1, F);
            times_pk1_mat(j) = toc(tic);
        
            % Compute PK2 Material response
            tic
            [stress, dtan ] = material_CSWP_PK2_hyperelasticity( dof, mat_pk2, F );
            times_pk2_mat(j) = toc(tic);
    
            % Compute Core Responses
            [EK_pk2, R_pk2, time_pk2] = core_pk2(k0, F, stress, dtan, nnElem, nn, N, ders, k0x, lambda, mu, M, e, Theta);
            [EK_pk1, R_pk1, time_pk1] = core_pk1(pk1, A, nnElem, nn, N, ders, k0x, lambda, mu, M, e, Theta);
            times_pk2_core(j) = time_pk2;
            times_pk1_core(j) = time_pk1;
        end
    
    
        % Append the Core Time data
        time_pk1_max_val(i) = max(times_pk1_core);
        time_pk1_min_val(i) = min(times_pk1_core);
        time_pk1_mean_val(i) = mean(times_pk1_core);
        time_pk1_std_val(i) = std(times_pk1_core);
        time_pk2_max_val(i) = max(times_pk2_core);
        time_pk2_min_val(i) = min(times_pk2_core);
        time_pk2_mean_val(i) = mean(times_pk2_core);
        time_pk2_std_val(i) = std(times_pk2_core);

        % Append the Mat Time data
        time_pk1_mat_max_val(i) = max(times_pk1_mat);
        time_pk1_mat_min_val(i) = min(times_pk1_mat);
        time_pk1_mat_mean_val(i) = mean(times_pk1_mat);
        time_pk1_mat_std_val(i) = std(times_pk1_mat);
        time_pk2_mat_max_val(i) = max(times_pk2_mat);
        time_pk2_mat_min_val(i) = min(times_pk2_mat);
        time_pk2_mat_mean_val(i) = mean(times_pk2_mat);
        time_pk2_mat_std_val(i) = std(times_pk2_mat);
    end

    % Visualize the Required Time Distributions
    n_elements = ipt_elevs;
    fig = figure("Color","w", "Position", [100, 100, 700, 600]);
    font_name = 'Helvetica';
    font_size = 11;
    
    subplot(2,1,1)
    hold on
    grid on
    set(gca, 'XScale', 'log', 'YScale', 'log', 'FontName', font_name, 'FontSize', font_size, 'LineWidth', 1.0);
    
    offset_factor = 1.03; 
    x_pk1 = n_elements / offset_factor;
    x_pk2 = n_elements * offset_factor;
    
    p1 = plot(n_elements, time_pk1_mean_val, '--k', 'LineWidth', 1.8);
    p2 = plot(n_elements, time_pk2_mean_val, '-b', 'LineWidth', 1.8);
    
    eb1_minmax = errorbar(x_pk1, time_pk1_mean_val, time_pk1_mean_val - time_pk1_min_val, time_pk1_max_val - time_pk1_mean_val, 'k', 'LineStyle', 'none', 'LineWidth', 1.2);
    eb2_minmax = errorbar(x_pk2, time_pk2_mean_val, time_pk2_mean_val - time_pk2_min_val, time_pk2_max_val - time_pk2_mean_val, 'b', 'LineStyle', 'none', 'LineWidth', 1.2);
    
    scatter(x_pk1, time_pk1_mean_val, 40, 'k', 'filled');
    scatter(x_pk2, time_pk2_mean_val, 40, 'b', 'filled');
    
    set(gca, 'XTick', n_elements, 'XTickLabel', string(n_elements));
    xlim([min(n_elements)/1.5, max(n_elements)*1.5]);
    
    global_min = min([time_pk1_min_val, time_pk2_min_val]);
    global_max = max([time_pk1_max_val, time_pk2_max_val]);
    ylim([global_min / 2, global_max * 2]);
    
    xlabel('Polynomial Degree', 'FontName', font_name, 'FontSize', font_size+1, 'FontWeight', 'bold');
    ylabel('Computation Time Core [s]', 'FontName', font_name, 'FontSize', font_size+1, 'FontWeight', 'bold');
    title('Computational Demand Comparison', 'FontName', font_name, 'FontSize', font_size+3, 'FontWeight', 'bold');
    
    legend([p1, p2], {'Arora-based PK1 formulation (Mean)', 'Novel PK2 formulation (Mean)'}, ...
        'Location', 'northeast', 'FontName', font_name, 'FontSize', font_size);
        
    hold off
end











if plot_K_computation_comparison
    geo_options.RefinementX = 2;
    geo_options.RefinementY = 2;
    plate = geo_square([0,0], 1, geo_options);
    mesh = build_iga_mesh(plate);
    u = zeros(mesh.nElems * 16 + 6, 1);
    eltype = 30;

    eps0 = [0.02, 0.03, 0.06]';
    k0 = [0.01,0.02,0.1]';
    

    % Store indicees for compared material models
    index_SVK_pk1 = 14; % Saint-Venant Kirchhoff with PK1
    index_SVK_pk2 = 114; % Saint-Venant Kirchhoff with PK2
    mat = default_mat();
    mat_pk1 = mat;
    mat_pk1.index = index_SVK_pk1;
    mat_pk2 = mat;
    mat_pk2.index = index_SVK_pk2;
    


    [K_pk1,r_pk1,times_pk1] = globalstiffness_CSWP_PK1_Arora_timed(eltype, plate, mesh, mat_pk1, u, 1, eps0, k0);
    [K_pk2,r_pk2,times_pk2] = globalstiffness_CSWP_PK2_timed(eltype, plate, mesh, mat_pk2, u, 1, eps0, k0);
    
    font_name = 'Helvetica';
    font_size = 11;
    
    figure('Color', 'w', 'Position', [100, 100, 600, 450]);
    subplot(1,3,1);
    hold on; grid on;
    set(gca, 'FontName', font_name, 'FontSize', font_size, 'LineWidth', 1.0);
    
    plot(times_pk1.K_times_core, 'r', 'LineWidth', 1.5);
    plot(times_pk2.K_times_core, 'b', 'LineWidth', 1.5);
    
    xlabel('Index', 'FontName', font_name, 'FontSize', font_size+1, 'FontWeight', 'bold');
    ylabel('Core Computation Time [s]', 'FontName', font_name, 'FontSize', font_size+1, 'FontWeight', 'bold');
    title('K{IJ} & R_{I} Computation Time Comparison', 'FontName', font_name, 'FontSize', font_size+2, 'FontWeight', 'bold');
    legend({'PK1 Formulation', 'PK2 Formulation'}, 'Location', 'northeast', 'FontName', font_name, 'FontSize', font_size);
    hold off;

    % Plot Mat computation time
    subplot(1,3,2)
    hold on; grid on;
    set(gca, 'FontName', font_name, 'FontSize', font_size, 'LineWidth', 1.0);
    
    plot(times_pk1.K_times_mat, 'r', 'LineWidth', 1.5);
    plot(times_pk2.K_times_mat, 'b', 'LineWidth', 1.5);
    
    xlabel('Index', 'FontName', font_name, 'FontSize', font_size+1, 'FontWeight', 'bold');
    ylabel('Constitutive Response Computation Time [s]', 'FontName', font_name, 'FontSize', font_size+1, 'FontWeight', 'bold');
    title('K{IJ} & R_{I} Computation Time Comparison', 'FontName', font_name, 'FontSize', font_size+2, 'FontWeight', 'bold');
    legend({'PK1 Formulation', 'PK2 Formulation'}, 'Location', 'northeast', 'FontName', font_name, 'FontSize', font_size);
    hold off;

end






%%%%%%%%%%%%%%%%%%%%%%%%%%


function [EK, RE, pk2_time] = core_pk2(k0, F, stress, dtan, nnElem, nn, N, ders, k0x, lambda, mu, M, e, Theta)
    tic;
    BN = zeros(6, nn*3);
    BG = zeros(6,nn,nn);

    for i = 1:nn
        % Corresponds to Eq. 52 in (1)
        BN(:,i*3-2:i*3) = [ F(1,1)*ders(1,i)     F(2,1)*ders(1,i)      F(3,1)*ders(1,i);
            F(1,2)*ders(2,i)     F(2,2)*ders(2,i)      F(3,2)*ders(2,i);
            N(i)*(k0(3)*F(2,3)-k0(2)*F(3,3))    N(i)*(k0(1)*F(3,3)-k0(3)*F(1,3))       N(i)*(k0(2)*F(1,3)-k0(1)*F(2,3)) ;
            F(1,1)*ders(2,i)+ F(1,2)*ders(1,i)  F(2,1)*ders(2,i)+F(2,2)*ders(1,i)   F(3,1)*ders(2,i)+F(3,2)*ders(1,i);
            (F(1,3)*ders(2,i) + N(i)*(k0(3)*F(2,2)-k0(2)*F(3,2)))  (F(2,3)*ders(2,i) + N(i)*(k0(1)*F(3,2)-k0(3)*F(1,2)))   (F(3,3)*ders(2,i) + N(i)*(k0(2)*F(1,2)-k0(1)*F(2,2)));
            (F(1,3)*ders(1,i) + N(i)*(k0(3)*F(2,1)-k0(2)*F(3,1)))  (F(2,3)*ders(1,i) + N(i)*(k0(1)*F(3,1)-k0(3)*F(1,1)))   (F(3,3)*ders(1,i) + N(i)*(k0(2)*F(1,1)-k0(1)*F(2,1))) ];

        for j=1:nn
            % Corresponds to Eq. 70 in (1)
            % Attention: kappa-related terms are reintroduced later in term G
            BG(:,i,j) = [ders(1,i)*ders(1,j);
                ders(2,i)*ders(2,j);
                -N(i)*N(j);                                
                ders(1,i)*ders(2,j)+ders(2,i)*ders(1,j);
                ders(2,i)*N(j)-N(i)*ders(2,j);            
                ders(1,i)*N(j)-N(i)*ders(1,j)   
                ];
        end
        RE2(3*i-2:3*i) = N(i)*(lambda+M*mu);
        EB(:,3*i-2:3*i) = N(i)*e;
        EC(:,3*i-2:3*i) = N(i)*M';
    end
    
    G = zeros(nnElem);
    % Corresponds to Eq. (71)
    for ii=1:nn
        for jj=1:nn
            G(ii*3-2:ii*3,jj*3-2:jj*3) = eye(3)*(stress(1)*BG(1,ii,jj)+stress(2)*BG(2,ii,jj)+stress(4)*BG(4,ii,jj))...
                +stress(3)*BG(3,ii,jj)*(k0x*k0x)+(stress(5)*BG(5,ii,jj)+stress(6)*BG(6,ii,jj))*k0x;
        end
    end
    
    RE = BN'*stress+RE2';
    EK = BN'*dtan*BN + G + kron(N(:)*N(:)',Theta);
    pk2_time = toc(tic);
end


function [EK, RE, pk1_time] = core_pk1(pk1, A, nnElem, nn, N, ders, k0x, lambda, mu, M, e, Theta)
    % ------------------------------------------------------------
    % Precompute tangent slices once per Gauss point
    % ------------------------------------------------------------
    tic
    A_i3_j3 = squeeze(A(:,3,:,3));

    A_i3_j1 = squeeze(A(:,3,:,1));
    A_i3_j2 = squeeze(A(:,3,:,2));

    A_i1_j3 = squeeze(A(:,1,:,3));
    A_i2_j3 = squeeze(A(:,2,:,3));

    A_i1_j1 = squeeze(A(:,1,:,1));
    A_i1_j2 = squeeze(A(:,1,:,2));
    A_i2_j1 = squeeze(A(:,2,:,1));
    A_i2_j2 = squeeze(A(:,2,:,2));

    % PK1 stress vector, column-wise:
    % Pvec = [pk1(:,1); pk1(:,2); pk1(:,3)]
    Pvec = pk1(:);

    RE  = zeros(nnElem, 1);
    RE2 = zeros(nnElem, 1);
    EK  = zeros(nnElem, nnElem);

    EB = zeros(3, nnElem);
    EC = zeros(3, nnElem);

    % ------------------------------------------------------------
    % Residual: Arora Eq. (39)
    %
    % H_J = N_J (P e3 x k0) + P e_alpha N_{J,alpha}
    %
    % Since P e3 x k0 = pk1(:,3) x k0 = -k0x * pk1(:,3)
    % ------------------------------------------------------------
    for J = 1:nn
        rows = 3*(J-1) + (1:3);
        NJ  = N(J);
        NJ1 = ders(1,J);
        NJ2 = ders(2,J);

        QJ = [ ...
            NJ1 * e, ...
            NJ2 * e, ...
           -NJ  * k0x ];

        RE(rows) = QJ * Pvec;

        % Constraint residual and Lagrange multiplier couplings
        RE2(rows) = NJ * (lambda + M * mu);

        EB(:,rows) = NJ * e;
        EC(:,rows) = NJ * M.';
    end

    % ------------------------------------------------------------
    % Stiffness: optimized Arora A_JI block assembly, Eq. (37)
    % ------------------------------------------------------------
    for J = 1:nn
        rows = 3*(J-1) + (1:3);
        NJ  = N(J);
        NJ1 = ders(1,J);
        NJ2 = ders(2,J);
    
        for I = 1:nn
            cols = 3*(I-1) + (1:3);
            NI  = N(I);
            NI1 = ders(1,I);
            NI2 = ders(2,I);
    
            % ----------------------------------------------------
            % Combined terms T1 + T2:
            %
            % T1 = -k0x * A_i3_j3 * k0x * NI * NJ
            % T2 = -k0x * (A_i3_j1*NI1 + A_i3_j2*NI2) * NJ
            %
            % Common factor: -k0x * (...) * NJ
            % ----------------------------------------------------
            T12 = -k0x * ( ...
                      A_i3_j3 * k0x * NI ...
                    + A_i3_j1 * NI1 ...
                    + A_i3_j2 * NI2 ) * NJ;
    
            % ----------------------------------------------------
            % Combined terms T3 + T4:
            %
            % T3 = (NJ1*A_i1_j3 + NJ2*A_i2_j3) * k0x * NI
            %
            % T4 = NJ1*(A_i1_j1*NI1 + A_i1_j2*NI2)
            %    + NJ2*(A_i2_j1*NI1 + A_i2_j2*NI2)
            %
            % Common factors: NJ1 and NJ2
            % ----------------------------------------------------
            T34 = NJ1 * ( ...
                      A_i1_j3 * k0x * NI ...
                    + A_i1_j1 * NI1 ...
                    + A_i1_j2 * NI2 ) ...
                + NJ2 * ( ...
                      A_i2_j3 * k0x * NI ...
                    + A_i2_j1 * NI1 ...
                    + A_i2_j2 * NI2 );
    
            EK(rows,cols) = T12 + T34;
        end
    end

    % Constraint stiffness:
    % block (J,I) = N(J) * N(I) * Theta
    EK = EK + kron(N(:) * N(:).', Theta);

    % Assemble residuals
    RE = (RE + RE2);
    pk1_time = toc(tic);
end