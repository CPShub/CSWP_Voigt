function plot_u_over_time(mesh, geo, us, title_text, options, yd)
    % 0. Standard-Optionen setzen
    if isempty(options)
        options = struct('grid', 1, 'yielding', 0, 'arrows', 1);
    end
    
    has_yd = 0;%(nargin >= 5 && ~isempty(yd));
    us = squeeze(us);
    [num_dofs, num_timesteps] = size(us);
    x0 = mesh.coords(:, 1); y0 = mesh.coords(:, 2); z0 = mesh.coords(:, 3);
    
    % UI Figure Setup
    fig = uifigure('Name', title_text, 'Color', 'w');
    ax = uiaxes(fig, 'Position', [70 220 850 520]);
    
    % State Management
    state = struct('q_ptr', [], 'p_def', [], 'g_ptr', [], ...
                   'yd_ptr_neg', [], 'yd_ptr_pos', [], 'p_ref', [], ...
                   'current_t', 1, 'current_scale', 1, 'opts', options);
    
    hold(ax, 'on'); grid(ax, 'off'); axis(ax, 'equal'); view(ax, 3);
    
    % --- CONTROL BUTTONS (TOP RIGHT) ---
    btn_grid = uibutton(fig, 'state', 'Text', 'Grid', 'Position', [70 140 80 25], ...
        'Value', state.opts.grid, 'ValueChangedFcn', @(b,e) toggleOpt('grid', b.Value));
    
    btn_arrows = uibutton(fig, 'state', 'Text', 'Arrows', 'Position', [160 140 80 25], ...
        'Value', state.opts.arrows, 'ValueChangedFcn', @(b,e) toggleOpt('arrows', b.Value));
    
    yd_enable_state = 'off';
    if has_yd
        yd_enable_state = 'on';
        u0 = zeros(size(us, 1), 1);
        yd_initial = get_integration_positions(geo, mesh, u0);
        yd_initial = [yd_initial, zeros(size(yd_initial, 1), 1)];
    end
    btn_yd = uibutton(fig, 'state', 'Text', 'Yielding', 'Position', [250 140 80 25], ...
        'Value', state.opts.yielding, 'Enable', yd_enable_state, ...
        'ValueChangedFcn', @(b,e) toggleOpt('yielding', b.Value));
    % --- TIME SLIDER & BUTTONS ---
    uilabel(fig, 'Position', [450 115 100 20], 'Text', 'Time Step', 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    
    slider_enable = 'on';
    slider_limits = [1, num_timesteps];
    if num_timesteps == 1
        slider_enable = 'off';
        slider_limits = [1, 1.001];
    end
    
    uibutton(fig, 'Text', '<', 'Position', [110 90 30 22], 'Enable', slider_enable, 'ButtonPushedFcn', @(~, ~) stepTime(-1));
    sld_time = uislider(fig, 'Position', [150 100 700 3], 'Limits', slider_limits, 'Value', 1, 'Enable', slider_enable, ...
        'ValueChangedFcn', @(s, e) updateTime(round(s.Value)));
    uibutton(fig, 'Text', '>', 'Position', [860 90 30 22], 'Enable', slider_enable, 'ButtonPushedFcn', @(~, ~) stepTime(1));
    % --- SCALE INPUT ---
    uilabel(fig, 'Position', [410 40 100 20], 'Text', 'Scale Factor:', 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
    ef_scale = uieditfield(fig, 'numeric', 'Position', [520 40 80 22], 'Value', 1, ...
        'ValueChangedFcn', @(ef, e) updateScale(ef.Value));
    renderFrame();
    % --- CALLBACKS ---
    function toggleOpt(field, val)
        state.opts.(field) = val;
        renderFrame();
    end
    function stepTime(direction)
        newVal = state.current_t + direction;
        if newVal >= 1 && newVal <= num_timesteps
            sld_time.Value = newVal;
            updateTime(newVal);
        end
    end
    function updateTime(t)
        state.current_t = t;
        renderFrame();
    end
    function updateScale(s)
        state.current_scale = s;
        renderFrame();
    end
    function renderFrame()
        % Daten extrahieren
        u_all = us(:, state.current_t);
        ux = u_all(1:3:end-6); uy = u_all(2:3:end-6); uz = u_all(3:3:end-6);
        x_def = x0 + ux * state.current_scale;
        y_def = y0 + uy * state.current_scale;
        z_def = z0 + uz * state.current_scale;
        
        % Aufräumen
        delete([state.q_ptr, state.p_def, state.g_ptr, state.yd_ptr_neg, state.yd_ptr_pos, state.p_ref]);
        
        % 1. Yielding Punkte
        if has_yd && state.opts.yielding
            % Determine the delta to the initial integration position
            d_yd = yd(:, :, state.current_t) - yd_initial;
            
            c_yd = yd_initial + d_yd .* state.current_scale;
            idx_pos = c_yd(:, 4) > 0;
            idx_neg = c_yd(:, 4) <= 0;
            state.yd_ptr_neg = scatter3(ax, c_yd(idx_neg,1), c_yd(idx_neg,2), c_yd(idx_neg,3), 10, [0.6 0.6 0.6], 'filled', 'MarkerFaceAlpha', 0.2);
            state.yd_ptr_pos = scatter3(ax, c_yd(idx_pos,1), c_yd(idx_pos,2), c_yd(idx_pos,3), 25, [1 0 0], 'filled');
        end
        % 2. Referenz & Pfeile
        if state.opts.arrows
            state.p_ref = scatter3(ax, x0, y0, z0, 15, [0.7 0.7 0.7], 'filled', 'MarkerFaceAlpha', 0.2);
            state.q_ptr = quiver3(ax, x0, y0, z0, ux*state.current_scale, uy*state.current_scale, uz*state.current_scale, 0, 'Color', [0.85 0.32 0.1], 'LineWidth', 0.8);
        end
        % 3. Grid
        if state.opts.grid
            n_pts = length(x_def); nz = 50; % Beispielwert, ggf. anpassen
            try
                X = reshape(x_def, [], nz); Y = reshape(y_def, [], nz); Z = reshape(z_def, [], nz);
                state.g_ptr = surf(ax, X, Y, Z, 'FaceColor', 'none', 'EdgeColor', [0.2 0.2 0.2], 'EdgeAlpha', 0.2, 'LineWidth', 3);
            catch
                tri = delaunay(x0, y0);
                state.g_ptr = trimesh(tri, x_def, y_def, z_def, 'Parent', ax, 'EdgeColor', [0.3 0.3 0.3], 'FaceColor', 'none', 'EdgeAlpha', 0.3, 'LineWidth', 3);
            end
        end
        % 4. Deformierte Punkte
        state.p_def = scatter3(ax, x_def, y_def, z_def, 100, [0 0.45 0.74], 'filled');
        
        title(ax, sprintf('%s | Step: %d | Scale: %.2f', title_text, state.current_t, state.current_scale));
    end
end