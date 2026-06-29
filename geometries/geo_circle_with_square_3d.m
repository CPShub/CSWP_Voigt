function all_nurbs = geo_circle_with_square_3d(center, radius, side_length, z_value)
if nargin < 3
    center = [0, 0];
    radius = 4;
    side_length = 2;
    z_value = 0; % Standardmäßig auf Z = 0
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

%  PATCH 1: Unten (Süden)
coefs1 = zeros(4,3,3,1); % Letzte Dimension = 1 für W-Richtung
coefs1(1:2, 1, 1, 1) = [-xc, -yc];          
coefs1(1:2, 2, 1, 1) = [0, -r_tangent];     
coefs1(1:2, 3, 1, 1) = [xc, -yc];
coefs1(1:2, 1, 2, 1) = [-mid_diag, -mid_diag]; 
coefs1(1:2, 2, 2, 1) = [0, -mid_r];         
coefs1(1:2, 3, 2, 1) = [mid_diag, -mid_diag];
coefs1(1:2, 1, 3, 1) = [-s, -s];            
coefs1(1:2, 2, 3, 1) = [0, -s];             
coefs1(1:2, 3, 3, 1) = [s, -s];
w_mat1 = ones(3,3); 
w_mat1(2,1) = w; 
all_nurbs{1} = finalize_patch_3d(coefs1, w_mat1, center, z_value);

%  PATCH 2: Rechts (Osten)
coefs2 = zeros(4,3,3,1);
coefs2(1:2, 1, 1, 1) = [s, -s];             
coefs2(1:2, 1, 2, 1) = [s, 0];              
coefs2(1:2, 1, 3, 1) = [s, s];
coefs2(1:2, 2, 1, 1) = [mid_diag, -mid_diag]; 
coefs2(1:2, 2, 2, 1) = [mid_r, 0];         
coefs2(1:2, 2, 3, 1) = [mid_diag, mid_diag];
coefs2(1:2, 3, 1, 1) = [xc, -yc];           
coefs2(1:2, 3, 2, 1) = [r_tangent, 0];      
coefs2(1:2, 3, 3, 1) = [xc, yc];
w_mat2 = ones(3,3); 
w_mat2(3,2) = w; 
all_nurbs{2} = finalize_patch_3d(coefs2, w_mat2, center, z_value);

%  PATCH 3: Oben (Norden)
coefs3 = zeros(4,3,3,1);
coefs3(1:2, 1, 1, 1) = [-s, s];             
coefs3(1:2, 2, 1, 1) = [0, s];              
coefs3(1:2, 3, 1, 1) = [s, s];
coefs3(1:2, 1, 2, 1) = [-mid_diag, mid_diag];  
coefs3(1:2, 2, 2, 1) = [0, mid_r];          
coefs3(1:2, 3, 2, 1) = [mid_diag, mid_diag];
coefs3(1:2, 1, 3, 1) = [-xc, yc];           
coefs3(1:2, 2, 3, 1) = [0, r_tangent];      
coefs3(1:2, 3, 3, 1) = [xc, yc];
w_mat3 = ones(3,3); 
w_mat3(2,3) = w; 
all_nurbs{3} = finalize_patch_3d(coefs3, w_mat3, center, z_value);

%  PATCH 4: Links (Westen)
coefs4 = zeros(4,3,3,1);
coefs4(1:2, 1, 1, 1) = [-xc, -yc];          
coefs4(1:2, 1, 2, 1) = [-r_tangent, 0];     
coefs4(1:2, 1, 3, 1) = [-xc, yc];
coefs4(1:2, 2, 1, 1) = [-mid_diag, -mid_diag]; 
coefs4(1:2, 2, 2, 1) = [-mid_r, 0];         
coefs4(1:2, 2, 3, 1) = [-mid_diag, mid_diag];
coefs4(1:2, 3, 1, 1) = [-s, -s];            
coefs4(1:2, 3, 2, 1) = [-s, 0];             
coefs4(1:2, 3, 3, 1) = [-s, s];
w_mat4 = ones(3,3); 
w_mat4(1,2) = w; 
all_nurbs{4} = finalize_patch_3d(coefs4, w_mat4, center, z_value);

%  PATCH 5: Zentrales Quadrat
c_sq = zeros(4,3,3,1);
x_vals = [-s, 0, s];
y_vals = [-s, 0, s];
for v = 1:3
    for u = 1:3
        c_sq(1:2, u, v, 1) = [x_vals(u), y_vals(v)] + center(:)';
        c_sq(3, u, v, 1)   = z_value; % Z-Koordinate setzen
        c_sq(4, u, v, 1)   = 1;
    end
end
% Erstellung des 3D-Patches mit Grad-0 Knotenvektor [0 1] für die W-Richtung
center_patch = nrbmak(c_sq, {[0 0 0 1 1 1], [0 0 0 1 1 1], [0 1]});
center_patch = nrbdegelev(center_patch, [1, 1, 0]); 
Ref = 3;
kv = 1/(Ref+1):1/(Ref+1):Ref/(Ref+1);
all_nurbs{5} = nrbkntins(center_patch, {kv, kv, []}); % Keine Knoten-Injektion in W-Richtung
end

function patch = finalize_patch_3d(coefs, weights, center, z_value)
    for u = 1:3
        for v = 1:3
            cw = weights(u,v);
            coefs(1:2, u, v, 1) = (coefs(1:2, u, v, 1) + center(:)) * cw;
            coefs(3, u, v, 1)   = z_value * cw; % Homogene Z-Koordinate gewichtet
            coefs(4, u, v, 1)   = cw;
        end
    end
    % Erstellung mit dem Grad 0 Knotenvektor [0 1] für die dritte Richtung
    patch = nrbmak(coefs, {[0 0 0 1 1 1], [0 0 0 1 1 1], [0 1]});
    patch = nrbdegelev(patch, [1, 1, 0]); 
    Ref = 3;
    kv = 1/(Ref+1):1/(Ref+1):Ref/(Ref+1);
    patch = nrbkntins(patch, {kv, kv, []}); % Keine Verfeinerung in W
end