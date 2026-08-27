function [F, F_c, J_c, J_e, I1] = chemo_kin(eps0, k0, x, dx_alpha, c, omega)
    % Compute relevant chemo-elastic kinematic quantities 
    e = eye(3);
    
    J_c = 1 + omega * c;
    F_c = J_c^(1/3) * e;

    F = (eps0 + cross(k0,x))*e(:,3)'+e+dx_alpha;
    J_e = det(F) / J_c;
    I1 = trace(det(F)^(-2/3) * (F'*F));
    
end
