function all_nurbs = geo_circle_with_square(center, radius, side_length, show_plot)
if nargin < 3
    center = [0, 0];
    radius = 4;
    side_length = 2;
    show_plot = 0;
end

s = side_length / 2; 
rad = pi/180;
w = cos(45*rad);
r_tangent = radius / w;
xc = radius * sin(45*rad);
yc = radius * cos(45*rad);
mid_r = (radius + s) / 2;
mid_diag = (xc + s) / 2;
degelev = [1,1];
all_nurbs = cell(1, 5); 

%  PATCH 1: Unten (Süden)
%  U -> läuft nach rechts (+X) | V -> läuft nach oben (+Y, von außen nach innen)
coefs1 = zeros(4,3,3);
% v = 1: Äußerer Kreisbogen unten
coefs1(1:2, 1, 1) = [-xc, -yc];          
coefs1(1:2, 2, 1) = [0, -r_tangent];     
coefs1(1:2, 3, 1) = [xc, -yc];
% v = 2: Mittelschicht
coefs1(1:2, 1, 2) = [-mid_diag, -mid_diag]; 
coefs1(1:2, 2, 2) = [0, -mid_r];         
coefs1(1:2, 3, 2) = [mid_diag, -mid_diag];
% v = 3: Innenkante zum Quadrat
coefs1(1:2, 1, 3) = [-s, -s];            
coefs1(1:2, 2, 3) = [0, -s];             
coefs1(1:2, 3, 3) = [s, -s];

w_mat1 = ones(3,3); 
w_mat1(2,1) = w; % Gewicht am Scheitelpunkt des Bogens (v=1, u=2)
all_nurbs{1} = finalize_patch(coefs1, w_mat1, center, degelev);


%  PATCH 2: Rechts (Osten)
% U -> läuft nach rechts (+X, von innen nach außen) | V -> läuft nach oben (+Y)
coefs2 = zeros(4,3,3);
% u = 1: Innenkante zum Quadrat (links)
coefs2(1:2, 1, 1) = [s, -s];             
coefs2(1:2, 1, 2) = [s, 0];              
coefs2(1:2, 1, 3) = [s, s];
% u = 2: Mittelschicht
coefs2(1:2, 2, 1) = [mid_diag, -mid_diag]; 
coefs2(1:2, 2, 2) = [mid_r, 0];         
coefs2(1:2, 2, 3) = [mid_diag, mid_diag];
% u = 3: Äußerer Kreisbogen rechts
coefs2(1:2, 3, 1) = [xc, -yc];           
coefs2(1:2, 3, 2) = [r_tangent, 0];      
coefs2(1:2, 3, 3) = [xc, yc];

w_mat2 = ones(3,3); 
w_mat2(3,2) = w; % Bogen liegt bei u=3, Scheitelpunkt bei v=2
all_nurbs{2} = finalize_patch(coefs2, w_mat2, center, degelev);


%  PATCH 3: Oben (Norden)
%  U -> läuft nach rechts (+X) | V -> läuft nach oben (+Y, von innen nach außen)
coefs3 = zeros(4,3,3);
% v = 1: Innenkante zum Quadrat
coefs3(1:2, 1, 1) = [-s, s];             
coefs3(1:2, 2, 1) = [0, s];              
coefs3(1:2, 3, 1) = [s, s];
% v = 2: Mittelschicht
coefs3(1:2, 1, 2) = [-mid_diag, mid_diag];  
coefs3(1:2, 2, 2) = [0, mid_r];          
coefs3(1:2, 3, 2) = [mid_diag, mid_diag];
% v = 3: Äußerer Kreisbogen oben
coefs3(1:2, 1, 3) = [-xc, yc];           
coefs3(1:2, 2, 3) = [0, r_tangent];      
coefs3(1:2, 3, 3) = [xc, yc];

w_mat3 = ones(3,3); 
w_mat3(2,3) = w; % Bogen liegt bei v=3, Scheitelpunkt bei u=2
all_nurbs{3} = finalize_patch(coefs3, w_mat3, center, degelev);


%  PATCH 4: Links (Westen)
%  U -> läuft nach rechts (+X, von außen nach innen) | V -> läuft nach oben (+Y)
coefs4 = zeros(4,3,3);
% u = 1: Äußerer Kreisbogen links
coefs4(1:2, 1, 1) = [-xc, -yc];          
coefs4(1:2, 1, 2) = [-r_tangent, 0];     
coefs4(1:2, 1, 3) = [-xc, yc];
% u = 2: Mittelschicht
coefs4(1:2, 2, 1) = [-mid_diag, -mid_diag]; 
coefs4(1:2, 2, 2) = [-mid_r, 0];         
coefs4(1:2, 2, 3) = [-mid_diag, mid_diag];
% u = 3: Innenkante zum Quadrat (rechts)
coefs4(1:2, 3, 1) = [-s, -s];            
coefs4(1:2, 3, 2) = [-s, 0];             
coefs4(1:2, 3, 3) = [-s, s];

w_mat4 = ones(3,3); 
w_mat4(1,2) = w; % Bogen liegt bei u=1, Scheitelpunkt bei v=2
all_nurbs{4} = finalize_patch(coefs4, w_mat4, center, degelev);


%  PATCH 5: Zentrales Quadrat
%  U -> läuft nach rechts (+X) | V -> läuft nach oben (+Y)
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

Ref = 3;
kv = 1/(Ref+1):1/(Ref+1):Ref/(Ref+1);
all_nurbs{5} = nrbkntins(center_patch, {kv, kv});


%  Visualisierung
show_plot = 0;
if show_plot
    figure; 
    hold on; 
    axis equal; 
    grid on;
    for i = 1:5
        plot_nurbs(all_nurbs{i}, 0, 1);
    end
    title('O-Mesh: Absolute KOS-Gleichheit aller Parameterachsen');
    view(2);
end
end

function patch = finalize_patch(coefs, weights, center, degelev)
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
    Ref = 3;
    kv = 1/(Ref+1):1/(Ref+1):Ref/(Ref+1);
    patch = nrbkntins(patch, {kv, kv});
end