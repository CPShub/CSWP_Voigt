function D = plastic_init_setting(mat, ngp, eltype)
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
% Initialize history variables and elastic stiffness matrix
% 
% input:
%       mat:    material definition
%       ngp:    number of integration points
%       eltype: element types,  
%           10 - plane strain element
%           20 - solid element
% 
% SIGMA : Stress for rate-form plasticity
%       : Left Cauchy-Green tensor XB for multiplicative plasticity
% 
% XQ    :   eltype == 10
%               1-4 = Back stress alpha, 5 = Effective plastic strain
%           eltype == 20
%               1-6 = Back stress alpha, 7 = Effective plastic strain
% D     : Elastic stiffness matrix
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%


E   = mat(2);
nu  = mat(3);
lam = E*nu/(1+nu)/(1-2*nu);
mu  = E/2/(1+nu);

if eltype == 10      % plane element
    global SIGMA XQ CPInv ZETA;

    D = [lam+2*mu    lam        0 ;
         lam         lam+2*mu   0 ; 
         0           0          mu;];
     
    SIGMA = zeros(4,ngp);
    XQ    = zeros(5,ngp);

elseif eltype == 30 % CSWP Planar Element with 3 DOF
    D = [lam+2*mu lam      lam      0  0  0;
         lam      lam+2*mu lam      0  0  0;
         lam      lam      lam+2*mu 0  0  0;
         0        0        0        mu 0  0;
         0        0        0        0  mu 0;
         0        0        0        0  0  mu];
    
    if use_vinayak
        % Code for Vinayak Style Plasticity
        global FP XI ZETA EE HH
        % plastic Def. Gradient, Kinematic Hardening Strain-type tensor,
        % Isotropic hardening strain-stype tensor

        % These update for each ngp point
        FP = zeros(3,3,ngp);
        XI = zeros(ngp, 1);
        ZETA = zeros(ngp, 1);

        % Define constant fourth-order macroscopic mdouli of elasticity tensor EE
        % Results in Voigt-Notation
        EE = compute_EE(lam, mu);

        % Define constant fourth-order Hill tensor
        % Results in voigt-Notation
        % Material properties taken from Vinayak 2023
        plast_mats = [
            0.672790077, 0.672790077, 0.666666667, ... 
            
            0.785052599, 0.785052599, 1.012246821];
        HH = compute_HH(plast_mats);
    else
        global SIGMA XQ CPInv ZETA;

        % Code for Herrnbock Style Plasticitiy
        SIGMA = zeros(6,ngp);
        XQ    = zeros(7,ngp);
        CPInv = zeros(3,3,ngp);
    
        % Set initial CPInv State (ones along diagonal entries)
        CPInv(1,1, :) = 1;
        CPInv(2,2, :) = 1;
        CPInv(3,3, :) = 1;
    
        % Set initial internal Hardening Variable State (zero)
        ZETA = zeros(ngp, 1);
    end



elseif eltype == 20  % solid element
    global SIGMA XQ CPInv ZETA;

    D = [lam+2*mu lam      lam      0  0  0;
         lam      lam+2*mu lam      0  0  0;
         lam      lam      lam+2*mu 0  0  0;
         0        0        0        mu 0  0;
         0        0        0        0  mu 0;
         0        0        0        0  0  mu];
    
    SIGMA = zeros(6,ngp);
    XQ    = zeros(7,ngp);
end

end