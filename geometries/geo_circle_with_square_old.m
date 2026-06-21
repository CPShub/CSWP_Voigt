function all_nurbs = geo_circle_with_square_old(center, radius, side_length, show_plot)
if nargin < 3
    center = [0, 0];
    radius = 4;
    side_length = 2;
    show_plot = 0;
end
% 
% if ~isempty(options)
%     if isfield(options, "show_plot"):
%         show_plot = options.show_plot;
%     else
%         show_plot = 1;
%     end
% else
%     show_plot = 1;
% end


s = side_length / 2; 
rad = pi/180;
w = cos(45*rad);
all_nurbs = cell(1, 5); 

%figure; hold on;

for i = 1:4
    % --- Basis points (unweighted Cartesian coordinates) ---
    % v=1: Square edge
    b(:,:,1) = [-s, -s;  0, -s;  s, -s]';
    % v=2: Mid layer
    mid_r = (radius + s) / 2;
    b(:,:,2) = [-mid_r, -mid_r; 0, -mid_r; mid_r, -mid_r]';
    % v=3: Circular arc
    xc = radius * sin(45*rad);
    yc = radius * cos(45*rad);
    % Point 2 is the tangent intersection at radius/w
    b(:,:,3) = [-xc, -yc; 0, -radius/w; xc, -yc]';

    % Corresponding weights for the layers
    weights = ones(3,3);
    weights(2,3) = w; % Only the middle point of the outer arc is rational

    % Rotation matrix
    phi = (i-1) * pi/2;
    rot_mat = [cos(phi), -sin(phi); sin(phi), cos(phi)];
    
    rotated_coefs = zeros(4,3,3);
    for u = 1:3
        for v = 1:3
            % 1. Get raw Cartesian point
            p_raw = b(:,u,v);
            
            % 2. Rotate and Translate in Cartesian space
            p_geom = rot_mat * p_raw + center(:);
            
            % 3. Convert to homogeneous coordinates (x*w, y*w, z*w, w)
            current_w = weights(u,v);
            rotated_coefs(1:2, u, v) = p_geom * current_w;
            rotated_coefs(3, u, v)   = 0;
            rotated_coefs(4, u, v)   = current_w;
        end
    end

    patch = nrbmak(rotated_coefs, {[0 0 0 1 1 1], [0 0 0 1 1 1]});
    patch = nrbdegelev(patch, [1, 1]);
    
    Ref = 3;
    kv = 1/(Ref+1):1/(Ref+1):Ref/(Ref+1);
    patch = nrbkntins(patch, {kv, kv});
    
    all_nurbs{i} = patch;
    %plot_nurbs(patch, 0, 1);
end

% --- Center Square (identical logic) ---
c_sq = zeros(4,2,2);
% Corner points
pts = {[-s,-s], [s,-s], [-s,s], [s,s]};
for j = 1:4
    [u,v] = ind2sub([2,2], j);
    p_final = pts{j}' + center(:);
    c_sq(1:2, u, v) = p_final;
    c_sq(4, u, v) = 1;
end

center_patch = nrbmak(c_sq, {[0 0 1 1], [0 0 1 1]});
center_patch = nrbdegelev(center_patch, [2, 2]);
kv_c = 1/(Ref+1):1/(Ref+1):Ref/(Ref+1);
all_nurbs{5} = nrbkntins(center_patch, {kv_c, kv_c});

if show_plot
    figure
    hold on; axis equal; view(2);
    plot_nurbs(all_nurbs{1}, 0, 1);
    plot_nurbs(all_nurbs{2}, 0, 1);
    plot_nurbs(all_nurbs{3}, 0, 1);
    plot_nurbs(all_nurbs{4}, 0, 1);
    plot_nurbs(all_nurbs{5}, 0, 1);
    
    axis([-inf,inf,-inf,inf]);
    view(2);
end






end