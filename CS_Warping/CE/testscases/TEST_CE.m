



% Setup the chemoelastic simulation
eps0 = [0, 0, 0]';
k0 = [0.1, 0, 0.2]';

% Define crossection and material
cs_size = 1; % Unit Square
square = geo_square([0, 0], cs_size, 0);
mesh = build_iga_mesh(square);
mat = Abdullah_mat();
mat.index = 310; % Chemoelastic NH
dt = 1;

% Constants for the nliga call
eltype = 30;    % element type: 10 - plane strain element, 20 - solid element, 30 - CSWP element
dbc =[];
tbc=[];

% setup previous solution vector
nCpts = mesh.nCpts;
c0 = 0.2 * ones(nCpts, 1);

% Compute the first CSWP 
fout = 69; 
%nliga_return1 = nliga_returns( eltype, square, mesh, mat, dbc, tbc, fout1 ,v01,k01);


% Run the simulation
nl_result_CE = nliga_returns_CE( eltype, square, mesh, mat, dbc, tbc, fout, eps0, k0, c0, dt);