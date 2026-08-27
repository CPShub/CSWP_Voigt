function all_nurbs = geo_circle_with_square(center, radius, side_length, varargin)
% This function generates a circular composit mesh with a square at the
% center. This formulation can resolve singularities based on geometry and
% is recommended instead of `geo_circle()`.
% Input:
    % center        - (1,2) vectoring containing center location
    % radius        - Radius of the circle
    % side_length   - Side-Length of the inner square
% Output:
    % all_nurbs     - geometry structure with defined nurb vectors
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

    if nargin < 3
        center = [0, 0];
        radius = 4;
        side_length = 2;
        show_plot_bool = 0;
    end
    
    % Handle additional geometry generation and visualization options
    if ~isempty(varargin)
        options = varargin{1};
    
        if isfield(options, "Refinement")
            Ref = options.Refinement;
        else
            Ref = 3;
        end
    
        if isfield(options, "degelev")
            degelev = options.degelev;
        else
            degelev = [1,1];
        end
    
        if isfield(options, "show_plot")
            show_plot_bool = options.show_plot;
        else
            show_plot_bool = 0;
        end
    else
        Ref = 3;
        degelev = [1,1];
        show_plot_bool = 0;
    end

    s = side_length / 2; 
    rad = pi/180;
    w = cos(45*rad);
    r_tangent = radius / w;
    xc = radius * sin(45*rad);
    yc = radius * cos(45*rad);
    mid_r = (radius + s) / 2;
    mid_diag = (xc + s) / 2;
    all_nurbs = cell(1, 5); 
    
    % PATCH 1: Bottom (South)
    % U -> runs to the right (+X) | V -> runs upward (+Y, from outside to inside)
    coefs1 = zeros(4,3,3);
    % v = 1: Outer circular arc
    coefs1(1:2, 1, 1) = [-xc, -yc];          
    coefs1(1:2, 2, 1) = [0, -r_tangent];     
    coefs1(1:2, 3, 1) = [xc, -yc];
    % v = 2: Intermediate layer
    coefs1(1:2, 1, 2) = [-mid_diag, -mid_diag]; 
    coefs1(1:2, 2, 2) = [0, -mid_r];         
    coefs1(1:2, 3, 2) = [mid_diag, -mid_diag];
    % v = 3: Inner edge of the square
    coefs1(1:2, 1, 3) = [-s, -s];            
    coefs1(1:2, 2, 3) = [0, -s];             
    coefs1(1:2, 3, 3) = [s, -s];
    
    w_mat1 = ones(3,3); 
    w_mat1(2,1) = w; % Weight at the apex of the arc (v=1, u=2)
    all_nurbs{1} = finalize_patch(coefs1, w_mat1, center, degelev, Ref);
    
    
    % PATCH 2: Right (East)
    % U -> runs to the right (+X, from inside to outside) | V -> runs upward (+Y)
    coefs2 = zeros(4,3,3);
    % u = 1: Inner edge of the square (left)
    coefs2(1:2, 1, 1) = [s, -s];             
    coefs2(1:2, 1, 2) = [s, 0];              
    coefs2(1:2, 1, 3) = [s, s];
    % u = 2: Intermediate layer
    coefs2(1:2, 2, 1) = [mid_diag, -mid_diag]; 
    coefs2(1:2, 2, 2) = [mid_r, 0];         
    coefs2(1:2, 2, 3) = [mid_diag, mid_diag];
    % u = 3: Outer circular arc on the right
    coefs2(1:2, 3, 1) = [xc, -yc];           
    coefs2(1:2, 3, 2) = [r_tangent, 0];      
    coefs2(1:2, 3, 3) = [xc, yc];
    
    w_mat2 = ones(3,3); 
    w_mat2(3,2) = w; % Weight at the apex of the arc (v=2, u=3)
    all_nurbs{2} = finalize_patch(coefs2, w_mat2, center, degelev, Ref);
    
    
    % PATCH 3: Top (North)
    % U -> runs to the right (+X) | V -> runs upward (+Y, from inside to outside)
    coefs3 = zeros(4,3,3);
    % v = 1: Inner edge of the square
    coefs3(1:2, 1, 1) = [-s, s];             
    coefs3(1:2, 2, 1) = [0, s];              
    coefs3(1:2, 3, 1) = [s, s];
    % v = 2: Intermediate layer
    coefs3(1:2, 1, 2) = [-mid_diag, mid_diag];  
    coefs3(1:2, 2, 2) = [0, mid_r];          
    coefs3(1:2, 3, 2) = [mid_diag, mid_diag];
    % v = 3: Outer circular arc on the top
    coefs3(1:2, 1, 3) = [-xc, yc];           
    coefs3(1:2, 2, 3) = [0, r_tangent];      
    coefs3(1:2, 3, 3) = [xc, yc];
    
    w_mat3 = ones(3,3); 
    w_mat3(2,3) = w; % Weight at the apex of the arc (v=3, u=2)
    all_nurbs{3} = finalize_patch(coefs3, w_mat3, center, degelev, Ref);
    
    
    % PATCH 4: Left (West)
    % U -> runs to the right (+X, from outside to inside) | V -> runs upward (+Y)
    coefs4 = zeros(4,3,3);
    % u = 1: Outer circular arc on the left
    coefs4(1:2, 1, 1) = [-xc, -yc];          
    coefs4(1:2, 1, 2) = [-r_tangent, 0];     
    coefs4(1:2, 1, 3) = [-xc, yc];
    % u = 2: Intermediate layer
    coefs4(1:2, 2, 1) = [-mid_diag, -mid_diag]; 
    coefs4(1:2, 2, 2) = [-mid_r, 0];         
    coefs4(1:2, 2, 3) = [-mid_diag, mid_diag];
    % u = 3: Inner edge of the square (right)
    coefs4(1:2, 3, 1) = [-s, -s];             
    coefs4(1:2, 3, 2) = [-s, 0];              
    coefs4(1:2, 3, 3) = [-s, s];
    
    w_mat4 = ones(3,3); 
    w_mat4(1,2) = w; % Weight at the apex of the arc (v=2, u=1)
    all_nurbs{4} = finalize_patch(coefs4, w_mat4, center, degelev, Ref);
    
    
    % PATCH 5: Central square
    % U -> runs to the right (+X) | V -> runs upward (+Y)
    c_sq = zeros(4,3,3);
    x_vals = [-s, 0, s];
    y_vals = [-s, 0, s];
    for v = 1:3
        for u = 1:3
            c_sq(1:2, u, v) = [x_vals(u), y_vals(v)] + center(:)';
            c_sq(3, u, v)   = 0;
            c_sq(4, u, v)   = 1;
        end
    end
    center_patch = nrbmak(c_sq, {[0 0 0 1 1 1], [0 0 0 1 1 1]});
    center_patch = nrbdegelev(center_patch, degelev);
    
    kv = 1/(Ref+1):1/(Ref+1):Ref/(Ref+1);
    all_nurbs{5} = nrbkntins(center_patch, {kv, kv});
    
    
    % Visualization
    if show_plot_bool
        figure; 
        hold on; 
        axis equal; 
        grid on;
        for i = 1:5
            plot_nurbs(all_nurbs{i}, 0, 1);
        end
        view(2);
    end
end

function patch = finalize_patch(coefs, weights, center, degelev, Ref)
    for u = 1:3
        for v = 1:3
            cw = weights(u,v);
            coefs(1:2, u, v) = (coefs(1:2, u, v) + center(:)) * cw;
            coefs(3, u, v)   = 0;
            coefs(4, u, v)   = cw;
        end
    end
    patch = nrbmak(coefs, {[0 0 0 1 1 1], [0 0 0 1 1 1]});
    patch = nrbdegelev(patch, degelev);
    kv = 1/(Ref+1):1/(Ref+1):Ref/(Ref+1);
    patch = nrbkntins(patch, {kv, kv});
end