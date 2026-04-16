function plot_color_flat_combined(flag, filenameA, filenameB, options)
% Plot stress, strain or displacement data from two files as a combined 
% surf over the respective cross-section. A customizable diagonal splits
% the view into one segment per file.
% Input:
	% flag              - Visualized data selection flag
	% filenameA         - Filename correspopnding to a valid .msh file
	% filenameB         - Filename correspopnding to a valid .msh file
	% options           - Struct containing further options for visualizations
    %   .show_title     - (Binary) Show plot title
    %   .show_coords    - (Binary) Show coordinate system in the corner
    %   .show_ticks     - (Binary) Show axes ticks
    %   .fontsize       - (Int) Fontsize for Different elements
    %   .A_text         - (Struct) Detailing Text written over data from A
    %   .B_text         - (Struct) Detailing Text written over data from B
    %   .diagonal       - (Struct) Detailing the positions spanning the
    %                       dividing diagonal line
    %   .cdata_mmax     - (Float) Upper Limit for colorbar
    %   .cdata_mmin     - (Float) Lower Limit for colorbar
    %   .cb_decimals    - (Binary) Number of displayed colorbar decimals
    %   .limit          - (Float) Limit data to range[-limit, limit]
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



% Handle call options
if isfield(options, "show_title")
    show_title = options.show_title;
else
    show_title = 1;
end
if isfield(options, "show_coords")
    show_coords = options.show_coords.flag;
    show_coords_center = options.show_coords.center; % Relative to the max Position in X and y
else
    show_coords = 1;
    show_coords_center = [-0.95, -0.95];
end

if isfield(options, "fontsize")
    txtsize = options.fontsize;
else
    txtsize = 20;
end

if isfield(options, "show_ticks")
    show_ticks = options.show_ticks;
else
    show_ticks = 1;
end

if isfield(options, "A_text")
    show_A_text = 1;
    A_text_pos = options.A_text.pos;
    A_text_text = options.A_text.text;
    A_text_font = options.A_text.font;
    A_text_color = options.A_text.color;
else
   show_A_text = 0;
end

if isfield(options, "B_text")
    show_B_text = 1;
    B_text_pos = options.B_text.pos;
    B_text_text = options.B_text.text;
    B_text_font = options.B_text.font;
    B_text_color = options.B_text.color;
else
   show_B_text = 0; 
end

if isfield(options, "cb_decimals")
    cb_decimals = options.cb_decimals;
else
    cb_decimals = 0;
end

if isfield(options, "cdata_mmin")
    cdata_mmin = options.cdata_mmin;
    cdata_mmax = options.cdata_mmax;
else
    cdata_mmin = 0;
    cdata_mmax = 0;
end

if isfield(options, "diagonal")
    diagonalA = options.diagonal.A;
    diagonalB = options.diagonal.B;
else
    diagonalA = [-2,2];
    diagonalB = [2,-2];
end

% Read and Divide the patches
vmeshA = read_visual_mesh( filenameA );
vmeshB = read_visual_mesh( filenameB );


% Trim each vmesh via a (diagonal) line
[trimAA, ~] = trim_patch_by_line(vmeshA, diagonalA, diagonalB);
[~, trimBB] = trim_patch_by_line(vmeshB, diagonalA, diagonalB);

% Dictionary Entries
cb_labels = {'$u_1$ in mm'; '$u_2$ in mm'; '$u_3$ in mm'; 
    '$|\vec{U}|$ in mm'; '$S_{11}\,in\,\frac{\text{kN}}{\text{mm}^2}$'; 
    '$S_{22}\,in\,\frac{kN}{mm^2}$'; '$S_{12}\,in\,\frac{kN}{mm^2}$'; 
    'von Mises in GPa'};
titles = {'u_1', 'u_2', 'u_3', 'u Magnitude', 
    'S_{11}', 'S_{22}', 'S_{12}', 'Mises'};
labeloffset = [-2, -2, -1.5, -1, -1, -1, -1, 1.5];

% Create the figure
f=figure;
f.Position = [1 1 650 500];
hold on;
if show_ticks
    xlabel('x');
    ylabel('y');
end
axis equal;
axis off;


% Plot only one half each (trimAA, trimBB)for vmeshA and vmeshB
trims = [trimAA, trimBB];
for j = 2:length(trims(1).vertices)
    for trimid = 1:4
        vmesh = trims(mod(trimid-1, 2)+1);
    
        % Draw the Undeformed Configuration (gray)
        if trimid <= 2
            face = cell2mat(vmesh.face(1));
            maxnum = max(max(face));
            vertices = cell2mat(vmesh.vertices(1));
            trivertex = vertices(1:maxnum,:);
            
            len = length(trivertex);
            p0 = fill(trivertex(:, 1), trivertex(:, 2), rand(len, 1));
            set(p0, 'Faces', face)
            
            % Update Face color to be light grey
            set(p0,'FaceColor', '#A19F99');
            set(p0,'EdgeColor','none');
            
            msize = get_visual_mesh_size(vmesh);
            axis(msize);
            pause(1.0);

            l = max(vertices);
        else
            % Draw the Deformed Configuration (colored) for each load step
    
            face = cell2mat(vmesh.face(j));
            maxnum = max(max(face));
            vertices = cell2mat(vmesh.vertices(j));
            trivertex = vertices(1:maxnum,:);
            displacement = cell2mat(vmesh.displacement(j));
            stress = cell2mat(vmesh.stress(j));
            % flag - color map: 
            %   1-U1, 2-U2, 3-U3, 
            %   4-U magnitude, 
            %   5-S11, 6-S22, 7-S12, 
            %   8-mises
            if flag == 1
                cdata = displacement(1:maxnum,1);
            elseif flag == 2
                cdata = displacement(1:maxnum,2);
            elseif flag == 3
                cdata = displacement(1:maxnum,3);
            elseif flag == 4
                cdata = sqrt(sum(displacement(1:maxnum, :).^2, 2));
            elseif flag == 5
                cdata = stress(1:maxnum,1);
            elseif flag == 6
                cdata = stress(1:maxnum,2);
            elseif flag == 7
                cdata = stress(1:maxnum,3);
            elseif flag == 8
                cdata = von_mises( stress(1:maxnum,:) );
            end
        
            % Add the visualized data
            p = patch(trivertex(:, 1), trivertex(:, 2), rand(length(trivertex), 1));
            set(p, 'Faces', face);
            set(p,'FaceColor','interp','FaceVertexCData',cdata);
            set(p,'EdgeColor','none');
        end
    end

    % add Title and colorbar
    [mmin, mmax] = bounds(cdata);
    if mmin == mmax
        mmin = -1;
        mmax = 1;
    end
    
    % Handel Singularity of the triangular elements (i.e. circle center)
    if cdata_mmin ~= 0 && cdata_mmax ~= 0
        mmin = cdata_mmin;
        mmax = cdata_mmax;
    end

    mmid = round((mmax - mmin) / 2 + mmin);

    cmp = colormap('jet');
    cb = colorbar('XTickLabel',{num2str(round(mmin, cb_decimals)), num2str(round(mmid, cb_decimals)), num2str(round(mmax, cb_decimals))}, ...
        'XTick', sort([mmin, mmid, mmax]));
    set(cb, 'FontSize', txtsize+4);
    set(cb, 'TickLabelInterpreter', "latex");
    set(gca,'CLim', [mmin mmax]);
    cb.Label.String = cb_labels(flag);
    cb.Label.VerticalAlignment = 'bottom';
    cb.Label.HorizontalAlignment = 'center';
    cb.Label.Interpreter = 'latex';
    cb.Label.FontSize = txtsize;
    cb.Label.Position(1) = cb.Label.Position(1) + labeloffset(flag);
    cb.FontSize = txtsize;
    if show_title
        if isfield(options, "given_title")
            title(options.given_title, 'interpreter','latex', 'FontSize', 20);
        else
            title("$$\bf{Step\," + num2str(j) + ":\," + titles(flag) + "}$$", 'interpreter','latex', 'FontSize', 20);
        end
    end
    drawnow; 

    if ~show_ticks
        set(gca,'xtick',[])
        set(gca,'ytick',[])
    end

    pause(0.5);
end
pause(2.0);



% Add a Diagonal Divide
xLimits = get(gca,'XLim');  % Get the range of the x axis
yLimits = get(gca,'YLim');  % Get the range of the y axis


% Determine intersection positions between limits and diagonal endpoints
lim_x_A = xLimits(1 + (diagonalA(1) > 0)); 
lim_y_A = yLimits(1 + (diagonalA(2) > 0));
lim_x_B = xLimits(1 + (diagonalB(1) > 0)); 
lim_y_B = yLimits(1 + (diagonalB(2) > 0));

t_A = min(abs([lim_x_A / diagonalA(1), lim_y_A / diagonalA(2)]));
t_B = min(abs([lim_x_B / diagonalB(1), lim_y_B / diagonalB(2)]));

P_A = t_A * diagonalA;
P_B = t_B * diagonalB;

plot([P_A(1), P_B(1)], [P_A(2), P_B(2)], "black", "LineWidth", 11);
scatter([P_A(1), P_B(1)], [P_A(2), P_B(2)], 90, "black", "filled");

% Add the coordinate system
if show_coords
    center = show_coords_center .* l;
    arrowlength = 0.15; % Relative to the Box size

    box_xsize = xLimits(2) - xLimits(1);
    box_ysize = yLimits(2) - yLimits(1);
    e3_offset = -0.08;
    quiver(center(1), center(2), arrowlength * box_xsize, 0, ...
        'linewidth',2.5, 'color', "k", "MaxHeadSize", 0.7) % x axis arrow
    quiver(center(1), center(2), 0, arrowlength * box_ysize, ...
        'linewidth',2.5, 'color', "k", "MaxHeadSize", 0.7) % y axis arrow
    scatter(center(1), center(2), 35, 'filled', "black")
    text(center(1) + arrowlength * box_xsize, center(2), "$\bf{E_1}$", ...
        "Interpreter","latex", "FontSize", txtsize)
    text(center(1), center(2) + arrowlength * box_ysize + 0.01, "$\bf{E_2}$", ...
        "Interpreter","latex", "FontSize", txtsize)
    text(center(1) + e3_offset, center(2) + e3_offset, "$\bf{E_3}$", ...
        "Interpreter","latex", "FontSize", txtsize)
end




% Write the associated Names
if show_A_text
    A_text_pos = A_text_pos .* l;
    text(A_text_pos(1), A_text_pos(2), A_text_text, 'FontSize', A_text_font, 'Color', A_text_color, 'FontWeight', 'bold', 'Interpreter', 'latex');
end
if show_B_text
    B_text_pos = B_text_pos .* l;
    text(B_text_pos(1), B_text_pos(2), B_text_text, 'FontSize', B_text_font, 'Color', B_text_color, 'FontWeight', 'bold', 'Interpreter', 'latex');
end

end

