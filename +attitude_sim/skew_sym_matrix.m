function [v_skew] = skew_sym_matrix(v)
% Skew-symmetric Matrix ===================================================
% Description: This function takes an input vector and converts it to the
% equivalent skew-symmetric (cross-product operator) matrix.
%
% Inputs:
%   v - The vector (3x1 or 1x3)
%
% Outputs:
%   v_skew - Skew-symmetric matrix (3x3)
%
% Created by: Cory Fraser - JAN 12, 2023
% Last Edits: Cory Fraser - JAN 12, 2023
% Copyright(c) 2023 by Cory Fraser
% =========================================================================
%% Organize the skew-symmetric matrix

v_skew = [  0       -v(3)   v(2)
            v(3)    0       -v(1)
            -v(2)   v(1)    0    ];

end