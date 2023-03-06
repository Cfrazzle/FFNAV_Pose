function [q_output] = rotmat2quat(Rarr)
%FROTMAT2QPARTS - quaternion parts from a frame rotation (body fixed) matrix
% original format used 
% a = scalar part, b,,d = vector part
%changed = q0, b,c,d = q1-3
%   Copyright 2017 The MathWorks, Inc.

dt = class(Rarr);
num = size(Rarr,3);
q0 = zeros(num,1, dt);
q1 = zeros(num,1, dt);
q2 = zeros(num,1, dt);
q3 = zeros(num,1, dt);

for ii=1:num
    R = Rarr(:,:,ii);
    [q1(ii), q2(ii), q3(ii), q0(ii)] = solveR2Q(R);
end

%q_output = [q0 q1 q2 q3];
q_output = [q1 q2 q3 q0];

end

function [q1,q2,q3,q0] = solveR2Q(R)
% The formula here can be found in "An Introduction to the Mathematics
% and Methods of Astrodynamics" by Richard H Battin. The formula is
% originally attributed to Stanley Shepperd's "Quaternion from Rotation
% Matrix" in Journal of Guidance and Control, 1978.

%Note the formulas in the Battin are for point rotation. This is a frame
%rotation conversion. The matrix element differences do not match those in
%Kuipers which is also frame rotation. Hence differences  are inverted
%which corresponds to transposing the matrix i.e. converting to a frame
%rotation.

tr = trace(R);

dd = [tr; diag(R)];
psquared = 1+2.*dd - trace(R);
[pmax,idx] = max(psquared);

switch (idx)
    case 1
        %angle is not near 180.
        pa = sqrt(pmax);
        q0 = 0.5*pa;
        invpa = 0.5./pa;
        q1 = invpa.*(R(2,3) - R(3,2));
        q2 = invpa.*(R(3,1) - R(1,3));
        q3 = invpa.*(R(1,2) - R(2,1));

    case 2
        %b is biggest
        pb = sqrt(pmax);
        q1 = 0.5*pb;
        invpb = 0.5./pb;

        q0 = invpb.*(R(2,3) - R(3,2));
        q2 = invpb.*(R(1,2) + R(2,1));
        q3 = invpb.*(R(3,1) + R(1,3));
    case 3
        %c is biggest
        pc = sqrt(pmax);
        q2 = 0.5*pc;
        invpc = 0.5./pc;

        q0 = invpc.*(R(3,1) - R(1,3));
        q1 = invpc.*(R(1,2) + R(2,1));
        q3 = invpc.*(R(2,3) + R(3,2));
    otherwise %4
        %d is biggest
        pd = sqrt(pmax);
        q3 = 0.5*pd;
        invpd = 0.5./pd;
        q0 = invpd.*(R(1,2) - R(2,1));
        q1 = invpd.*(R(3,1) + R(1,3));
        q2 = invpd.*(R(2,3) + R(3,2));
end

%Make first part of the quaternion positive
if q0 < 0
    q0 = -q0;
    q1 = -q1;
    q2 = -q2;
    q3 = -q3;
end

end
