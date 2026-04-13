function vmesh = plot_visual_mesh2d(flag, filename, z_options)
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
% Plot the results of visualized mesh model
%  Input:
%    flag - color map: 1-U1, 2-U2, 3-U magnitude, 4-S11, 5-S22, 6-S12, 7-mises
%    filename - filename for read
%  Output: 
%    vmesh - visualized mesh structure
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

fname = filename;
vmesh = read_visual_mesh( fname );               % read the saved .msh file

% % make movie
figure;
% videoName = [str_temp, 'Output\', filename, '.avi'];
% aviobj = VideoWriter(videoName,'Uncompressed AVI');
% aviobj.FrameRate = 30;
% open(aviobj);
% set(fig, 'outerposition',get(0,'screensize'));


% Draw the Mesh Grid
face = cell2mat(vmesh.face(1));
maxnum = max(max(face));
vertices = cell2mat(vmesh.vertices(1));
trivertex = vertices(1:maxnum,:);
displacement = cell2mat(vmesh.displacement(1));
stress = cell2mat(vmesh.stress(1));  
linmesh = cell2mat(vmesh.linmesh(1));

if nargin == 2
    % only flat 2d Case relevant 
    show_z = 0;
    z_scale = 1;
    limi = 1e30;
else
    % Accept 3d options
    if isfield(z_options, "show_z")
        show_z = z_options.show_z;
    else
        show_z = 0;
    end
    
    if isfield(z_options, "scale_z")
        z_scale = z_options.scale_z;
    else
        z_scale = 100;
    end
    
    if isfield(z_options, "limit")
        limi = z_options.limit;
    else
        limi = 1e30;
    end
end

if flag == 1
    cdata = displacement(1:maxnum,1);
elseif flag == 2
    cdata = displacement(1:maxnum,2);
elseif flag == 3
    cdata = sqrt(displacement(1:maxnum,1).^2 + displacement(1:maxnum,2).^2);
elseif flag == 4
    cdata = stress(1:maxnum,1);
elseif flag == 5
    cdata = stress(1:maxnum,2);
elseif flag == 6
    cdata = stress(1:maxnum,3);
elseif flag == 7
    cdata = von_mises( stress(1:maxnum,:) );
elseif flag == 8
    cdata = displacement(1:maxnum,3);
end

% Prepare 3d mesh
len = length(trivertex);

% Update z value according to cdata
if show_z
    cdata_scaled = cdata * z_scale;
    cdata_clamped = min(max(cdata_scaled,-limi),limi);
    trivertex(:, 3) = cdata_clamped;
else
    trivertex(:, 3) = 0;
end

if show_z
    p = fill3(trivertex(:, 1), trivertex(:, 2), trivertex(:, 3), rand(len, 1));
    set(p, 'Faces', face)
else
    p = patch('Faces',face, 'Vertices', trivertex);
end

% Update Face color according to cdata
set(p,'FaceColor','interp','FaceVertexCData',cdata);
set(p,'EdgeColor','none');

hold on; % plot knot curves
aa = zeros(size(linmesh,1),1);
for j = 1:size(linmesh,1)
    vv = vertices(linmesh(j,:),:);
    aa(j) = plot(vv(:,1), vv(:,2), 'k-');
end    
title(['Step ',num2str(1)]);
hcb = colorbar;
xlabel('x');
ylabel('y');
% set(gcf,'color','black');
axis equal;
msize = get_visual_mesh_size(vmesh);
axis(msize);
box on;
if flag == 1,  title(hcb,'U_x');
elseif flag == 2,  title(hcb,'U_y');
elseif flag == 3,  title(hcb,'U Magnitude');
elseif flag == 4,  title(hcb,'S_{11}');
elseif flag == 5,  title(hcb,'S_{22}');
elseif flag == 6,  title(hcb,'S_{12}');
elseif flag == 7,  title(hcb,'Mises');
end

pause(1.0);


for j = 2:length(vmesh.vertices)
    face = cell2mat(vmesh.face(j));
    maxnum = max(max(face));
    vertices = cell2mat(vmesh.vertices(j));
    trivertex = vertices(1:maxnum,:);
    displacement = cell2mat(vmesh.displacement(j));
    stress = cell2mat(vmesh.stress(j));    
    linmesh = cell2mat(vmesh.linmesh(j));
    if flag == 1
        cdata = displacement(1:maxnum,1);
    elseif flag == 2
        cdata = displacement(1:maxnum,2);
    elseif flag == 3
        cdata = sqrt(displacement(1:maxnum,1).^2 + displacement(1:maxnum,2).^2);
    elseif flag == 4
        cdata = stress(1:maxnum,1);
    elseif flag == 5
        cdata = stress(1:maxnum,2);
    elseif flag == 6
        cdata = stress(1:maxnum,3);
    elseif flag == 7
        cdata = von_mises( stress(1:maxnum,:) );
        % Limit stress data to 1e30
        limi = 1e30;
        cdata = cdata > limi;
    elseif flag == 8
        cdata = displacement(1:maxnum,3);
    end

    % Update trivertex z entries
    if show_z
        cdata_scaled = cdata * z_scale;
        cdata_clamped = min(max(cdata_scaled,-limi),limi);
        trivertex(:, 3) = cdata_clamped;
    else
        trivertex(:, 3) = 0;
    end

    set(p, 'Faces',face, 'Vertices', trivertex, 'FaceVertexCData',cdata);
    for k = 1:size(linmesh,1)
        vv = vertices(linmesh(k,:),:);
        set(aa(k),'xdata', vv(:,1), 'ydata', vv(:,2));
    end       
    title(['Step ',num2str(j)]);
    hcb = colorbar;
    if flag == 1,  title(hcb,'U_x');
    elseif flag == 2,  title(hcb,'U_y');
    elseif flag == 3,  title(hcb,'U Magnitude');
    elseif flag == 4,  title(hcb,'S_{11}');
    elseif flag == 5,  title(hcb,'S_{22}');
    elseif flag == 6,  title(hcb,'S_{12}');
    elseif flag == 7,  title(hcb,'Mises');
    end
    drawnow; 

    pause(0.5);
end
pause(2.0);

end

