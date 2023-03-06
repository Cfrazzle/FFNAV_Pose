function qm = quat_mult(q1,q2)
%function qm=quat_mult(q1,q2);
%
% This program multiplies two quaternion sets.
%
%  The inputs are:
%    q1 = first quaternion set (mx4)
%    q1 = [q1_vec q2_vec q3_vec q0_scalar];
%    q2 = second quaternion set (mx4)
%
%  The output is q1xq2
% John L. Crassidis 4/24/95
qm(:,1) =  q1(:,4).*q2(:,1) + q1(:,3).*q2(:,2) - q1(:,2).*q2(:,3) + q1(:,1).*q2(:,4);
qm(:,2) = -q1(:,3).*q2(:,1) + q1(:,4).*q2(:,2) + q1(:,1).*q2(:,3) + q1(:,2).*q2(:,4);
qm(:,3) =  q1(:,2).*q2(:,1) - q1(:,1).*q2(:,2) + q1(:,4).*q2(:,3) + q1(:,3).*q2(:,4);
qm(:,4) = -q1(:,1).*q2(:,1) - q1(:,2).*q2(:,2) - q1(:,3).*q2(:,3) + q1(:,4).*q2(:,4);

qm = my_quaternion.quat_normalize(qm);

end