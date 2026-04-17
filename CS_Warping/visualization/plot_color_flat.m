function vmesh = plot_color_flat(flag, filename, options)
% Plot stress, strain or displacement data as a surf over the cross-section
% Input:
	% flag              - Visualized data selection flag
	% filename          - Filename correspopnding to a valid .msh file
	% options           - Struct containing further options for visualizations
    %   .show_title     - (Binary) Show plot title
    %   .show_coords    - (Binary) Show coordinate system in the corner
    %   .show_ticks     - (Binary) Show axes ticks
    %   .show_cb_title  - (Binary) Show colorbar title
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



% read the saved .msh file
fname = filename;
vmesh = read_visual_mesh( fname );               
figure;

% Handle call options
if isfield(options, "show_title")
    show_title = options.show_title;
else
    show_title = 1;
end
if isfield(options, "show_coords")
    show_coords = options.show_coords;
else
    show_coords = 1;
end
if isfield(options, "show_ticks")
    show_ticks = options.show_ticks;
else
    show_ticks = 1;
end

if isfield(options, "show_cb_title")
    show_cb_title = options.show_cb_title;
else
    show_cb_title = 0;
end

if isfield(options, "limit")
    limi = options.limit;
else
    limi = 1e30;
end

% Text and label positioning via dictionary entries
cb_labels = {'$u_1\,in\,mm$'; '$u_2\,in\,mm$'; '$u_3\,in\,mm$'; 
    '$u\,Magnitude\,in\,mm$'; '$S_{11}\,in\,\frac{N}{mm^2}$'; 
    '$S_{22}\,in\,\frac{N}{mm^2}$'; '$S_{12}\,in\,\frac{N}{mm^2}$'; 
    '$Mises\,in\,\frac{N}{mm^2}$'};
titles = {'u_1', 'u_2', 'u_3', 'u Magnitude', 'S_{11}', ...
    'S_{22}', 'S_{12}', 'Mises'};
labeloffset = [-2, -2, -1, -1, -1, -1, -1, -1];

% Constant Textsize
txtsize = 16;

% Draw the Undeformed Configuration underneath
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

hold on;
if show_ticks
    xlabel('x');
    ylabel('y');
end
axis equal;
msize = get_visual_mesh_size(vmesh);
axis(msize);
box on;

pause(1.0);

% Draw the Deformed Configuration (colored) for each load step
for j = 2:length(vmesh.vertices)
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
    len = length(trivertex);
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

    % Limit the cdata to [-limi, limi] range
    cdata = min(max(cdata,-limi),limi);


    % Add the visualized data
    p = patch(trivertex(:, 1), trivertex(:, 2), rand(len, 1));
    set(p, 'Faces', face);
    set(p,'FaceColor','interp','FaceVertexCData',cdata);
    set(p,'EdgeColor','none');

    % Add the coordinate system
    if show_coords
        center = [0,0];
        arrowlength = 0.15; % Relative to the Box size
        xLimits = get(gca,'XLim');  % Get the range of the x axis
        yLimits = get(gca,'YLim');  % Get the range of the y axis
    
        box_xsize = xLimits(2) - xLimits(1);
        box_ysize = yLimits(2) - yLimits(1);
        e3_offset = -0.06;
        quiver(center(1), center(2), arrowlength * box_xsize, 0, ...
            'linewidth',2.5, 'color', "k", "MaxHeadSize", 0.7) % x axis arrow
        quiver(center(1), center(2), 0, arrowlength * box_ysize, ...
            'linewidth',2.5, 'color', "k", "MaxHeadSize", 0.7) % y axis arrow
        scatter(center(1), center(2), 35, 'filled', "black")
        text(center(1) + arrowlength * box_xsize, center(2), "$\bf{E_1}$", ...
            "Interpreter","latex", "FontSize", txtsize)
        text(center(1), center(2) + arrowlength * box_ysize, "$\bf{E_2}$", ...
            "Interpreter","latex", "FontSize", txtsize)
        text(center(1) + e3_offset, center(2) + e3_offset, "$\bf{E_3}$", ...
            "Interpreter","latex", "FontSize", txtsize)
    end

    % Configure the Colorbar
    [mmin, mmax] = bounds(cdata);
    if mmin == mmax
        mmin = -1;
        mmax = 1;
    end
    cmp = colormap('jet');
    cb = colorbar('XTickLabel',{num2str(mmin),num2str(mmax)}, ...
        'XTick', sort([mmin, mmax]));
    if show_cb_title
        set(cb, 'FontSize', txtsize-2);
        set(cb, 'TickLabelInterpreter', "latex");
        set(gca,'CLim', [mmin mmax]);
        cb.Label.String = cb_labels(flag);
        cb.Label.VerticalAlignment = 'bottom';
        cb.Label.HorizontalAlignment = 'center';
        cb.Label.Interpreter = 'latex';
        cb.Label.FontSize = 20;
        cb.Label.Position(1) = cb.Label.Position(1) + labeloffset(flag);
    end

    % Set the title
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

end

