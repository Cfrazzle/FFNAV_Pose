function qm = quat_err(q1,q2)
%function qm=quat_err(q1,q2);
%
% This program determines the error quaternion q1*q2^-1. (q_actual*qcommand^-1)
%
%  The inputs are:
%    q1 = first quaternion set (mx4)
%    q2 = second quaternion set (mx4)
%
%  The output is dq = [dq_13_vec; dq4_scalar] = q1 x q2^-1
%  dq_13_vec = Theta^T(qc) * q_actual
%  dq_4_scalar = q_actual^T * q_command
%
% Reference: Crassidis, eq 7.2, page 289-290
% John L. Crassidis 4/24/95
%% Formula for Reference
% delta_q = [delta_qvec delta_q_scalar]
% with delta_q_vec = Sigma^T(qc) *qvec
% with delta_q_scalar = q^T * qc
% and Sigma = [q_scalar -qvec(3) qvec(2)...
%              qvec(3) q_Scalar -qvec(1)...
%             -qvec(2) qvec(1)q_scalar
%             -qvec(1) -qvec(2) -qvec(3) ]
%
% and Sigma^T = [ q_scalar qvec(3) -qvec(2) -qvec(1)...
%                -qvec(3) q_Scalar qvec(1) -qvec(2)
%                 qvec(2) -qvec(2) q_scalar -qvec(3)
%
% Other Implementation for a single set of quaternions
% quat_mult_mat = [qd(4) qd(3) -qd(2) -qd(1)
%     -qd(3) qd(4) qd(1) -qd(2)
%     qd(2) -qd(1) qd(4) -qd(3)
%     qd(1) qd(2) qd(3) qd(4)];
% 
% quat_err = quat_mult_mat * quat;

% Create the conjugate
q2_conj = q2;
q2_conj(:,1:3) = -q2(:,1:3);

% quatmultiply(Qact,quatinv(Qc))
qm = my_quaternion.quat_mult(q1,q2_conj);

%% Compare results using quaternion objects
Qact    = normalize(quaternion([q1(:,4) q1(:,1:3)]));
Qc      = normalize(quaternion([q2(:,4) q2(:,1:3)]));
Qerr    = normalize(quatmultiply(Qact,quatinv(Qc)));
Qm      = normalize(quaternion([qm(:,4) qm(:,1:3)]));

if ~all((norm(Qerr - Qm) < 1e-5))
    error('Quaternions don''t match')
end

end