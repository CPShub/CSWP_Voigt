function plot_u_quiver(mesh, u, nit)
    % 1. Extraction
    u_pos = u(1:end-6);
    x = mesh.coords(:, 1);
    y = mesh.coords(:, 2);
    z = mesh.coords(:, 3);
    
    u_pos_x = u_pos(1:3:end); u_pos_x = u_pos_x(:);
    u_pos_y = u_pos(2:3:end); u_pos_y = u_pos_y(:);
    u_pos_z = u_pos(3:3:end); u_pos_z = u_pos_z(:);
    
    % Calculate the magnitude of displacement
    u_mag = sqrt(u_pos_x.^2 + u_pos_y.^2 + u_pos_z.^2);
    [maxVal, maxIdx] = max(u_mag);
    [minVal, minIdx] = min(u_mag);

    if maxVal == minVal
        maxVal = 1;
        minVal = -1;
    end

    % Pre-calculate the surface interpolation
    [xq, yq] = meshgrid(linspace(min(x), max(x), 50), linspace(min(y), max(y), 50));
    zq = griddata(x, y, z, xq, yq, 'v4');
    uq = griddata(x, y, u_mag, xq, yq, 'v4');

    % Create UI Figure and Axes
    fig = uifigure('Name', sprintf('Iteration %d', nit), 'Position', [100 100 850 650]);
    ax = uiaxes(fig, 'Position', [50 120 750 500]);
    
    % Initial Plot Setup
    hold(ax, 'on');
    grid(ax, 'on');
    axis(ax, 'equal');
    view(ax, 3);
    %ax.XLim = [-0.75 0.75];
    %ax.YLim = [-0.75 0.75];
    %ax.ZLim = [-0.75 0.75];
    xlabel(ax, "X"); ylabel(ax, "Y"); zlabel(ax, "Z");
    title(ax, sprintf('Iteration %d', nit));

    % Color setup
    colormap(ax, jet);
    clim(ax, [minVal maxVal]); % Modern version of set(gca, 'CLim'...)
    
    % Map Max and Min values to the 'jet' colormap for the markers
    cmap = jet(256);
    % Find color index for Max (top of map) and Min (bottom of map)
    colorMax = cmap(end, :); 
    colorMin = cmap(1, :);

    % Plot static elements
    scatter3(ax, x, y, z, 10, 'blue', 'filled', 'DisplayName', 'Mesh Nodes');
    surf(ax, xq, yq, zq, uq, 'EdgeColor', 'none', 'FaceAlpha', 0.6, 'HandleVisibility', 'off');
    
    % Plot Max/Min Highlight Points with corresponding colormap colors
    plot3(ax, x(maxIdx), y(maxIdx), z(maxIdx), 'o', 'MarkerFaceColor', colorMax, 'MarkerEdgeColor', 'k', 'MarkerSize', 10, ...
        'DisplayName', sprintf('Max: %.2e', maxVal));
    plot3(ax, x(minIdx), y(minIdx), z(minIdx), 'o', 'MarkerFaceColor', colorMin, 'MarkerEdgeColor', 'k', 'MarkerSize', 10, ...
        'DisplayName', sprintf('Min: %.2e', minVal));

    % Place Legend at Top Left
    lgd = legend(ax, 'Location', 'northwest');

    cb = colorbar(ax);
    cb.Label.String = 'Displacement Magnitude |u|';

    % Create the Adaptive Slider
    sld = uislider(fig, ...
        'Position', [125 60 600 3], ...
        'Limits', [1 1000], ...
        'Value', 1, ...
        'ValueChangedFcn', @(sld, event) updatePlot(sld.Value));
    
    uilabel(fig, 'Position', [375 75 100 20], 'Text', 'Scale Factor', 'FontWeight', 'bold');

    % Initial Quiver Plot
    q_ptr = quiver3(ax, x, y, z, u_pos_x, u_pos_y, u_pos_z, 0, 'r', 'HandleVisibility', 'off');

    % Nested Update Function
    function updatePlot(currentScale)
        if ishandle(q_ptr)
            delete(q_ptr);
        end
        q_ptr = quiver3(ax, x, y, z, ...
                        u_pos_x * currentScale, ...
                        u_pos_y * currentScale, ...
                        u_pos_z * currentScale, 0, 'r', 'HandleVisibility', 'off');
    end
end
