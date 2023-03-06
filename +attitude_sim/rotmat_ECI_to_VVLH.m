function [C_LI] = rotmat_ECI_to_VVLH(r_sat, v_sat)

% From Crassidis
%     O_3I = -r_sat/norm(r_sat);
%     O_2I = -(cross(r_sat,v_sat))/norm((cross(r_sat,v_sat)));
%     O_1I = cross(O_2I,O_3I);
%     A_IO = [O_1I O_2I O_3I];

%Determine unit vectors for VVLH Frame
Lz      = -r_sat/norm(r_sat);
Ly      = -cross(r_sat, v_sat)/norm(cross(r_sat, v_sat));
Lx      = cross(Ly,Lz);

%Align unit vectors into rows
Lx      = Lx';
Ly      = Ly';
Lz      = Lz';

%Arranging vectrix for LVLH Frame
F_VVLH  = [ Lx
    Ly
    Lz ];

%ECI Frame is defined as such
F_ECI   = [ 1 0 0
    0 1 0
    0 0 1 ];

%Compute roration matrix - ECI to LVLH
C_LI    = F_VVLH*F_ECI';

end