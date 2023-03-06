 function [r_rel_LVLH, v_rel_LVLH] = To_LVLH(rt, r_rel, vt, v_rel)
% =========================================================================
% Description: This function converts the relative r and v from ECI Frame
% to LVLH Frame by applying the appropriate rotation matrices.
% =========================================================================
%Determine unit vectors for LVLH Frame
Lx      = rt/norm(rt);
Lz      = cross(rt, vt)/norm(cross(rt, vt));
Ly      = cross(Lz, Lx);

%Align unit vectors into rows
Lx      = Lx';
Ly      = Ly';
Lz      = Lz';

%Arranging vectrix for LVLH Frame
F_LVLH  = [ Lx
            Ly
            Lz ];

%ECI Frame is defined as such
F_ECI   = [ 1 0 0 
            0 1 0
            0 0 1 ];
     
%Compute roration matrix - ECI to LVLH
C_LI    = F_LVLH*F_ECI';

%Computing relative velocity in LVLH frame, accounting for rotation 
omega =  cross(rt,vt)/norm(rt)^2;
v_rel = v_rel - cross(omega,r_rel);

%Apply rotation to obtain relative LVLH parameters
r_rel_LVLH = C_LI*r_rel; 
v_rel_LVLH = C_LI*v_rel;

end
% =========================================================================