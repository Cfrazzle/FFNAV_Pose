function qo = mySlerp(q1, q2, t)
%  QO = SLERP(Q1,Q2,T) spherically interpolates between Q1 and Q2 by the
%  interpolation coefficient T. T is a single or double precision number
%  between 0 and 1, inclusive.  The inputs Q1, Q2, and T must have
%  compatible sizes. In the simplest cases, they can be the same size or
%  any one can be a scalar.  Two inputs have compatible sizes if, for every
%  dimension, the dimension sizes of the inputs are either the same or one
%  of them is 1.
%
%  % Example:
%       e = deg2rad([40 20 10; 50 10 5; 45 70 1]);
%       q = quaternion(e, 'euler', 'ZYX', 'frame');
%       qs = slerp(q(1), q(2), 0.7);
%       rad2deg(euler(qs, 'ZYX', 'frame'))
%
%   See also QUATERNION, MEANROT
%   Copyright 2018-2021 The MathWorks, Inc.
%q1 and q2 and t must be compatible sizes and vectors
%% Algorithm
% Normalize an expand
% q1n = normalize(q1) .* ones(size(q2), 'like', q2);
% q2n = normalize(q2) .* ones(size(q1), 'like', q1);
q1n = my_quaternion.quat_normalize(q1) .* ones(size(q2), 'like', q2);
q2n = my_quaternion.quat_normalize(q2) .* ones(size(q1), 'like', q1);

% Get parts
%[a1, b1, c1, d1] = parts(q1n);
%[a2, b2, c2, d2] = parts(q2n);
[a1, b1, c1, d1] = deal(q1n(:,1), q1n(:,2), q1n(:,3), q1n(:,4));
[a2, b2, c2, d2] = deal(q2n(:,1), q2n(:,2), q2n(:,3), q2n(:,4));

% Implement quaternion dot product, inline
dp = a1.*a2 + b1.*b2 + c1.*c2 + d1.*d2;

% Negative dot product, the quaternions aren't pointing the same way (one
% pos, one negative). Flip the second one.
if isempty(coder.target)
    dpidx = dp < 0;
    if any(dpidx(:))
        q2n(dpidx,:) = -q2n(dpidx,:);
        dp(dpidx) = -dp(dpidx);
    end
    dp(dp > 1) = 1;
else
    for ii=1:numel(dp)
        if dp(ii) < 0
            q2n(ii,:) = -q2n(ii,:);
            dp(ii) = -dp(ii);
        end
        if dp(ii) > 1
            dp(ii) = 1;
        end
    end
end
theta0 = acos(dp);

sinv = 1./sin(theta0);
qnumerator = q1n.*sin((1- t).*theta0) + q2n.*sin(t.*theta0);
qo =  qnumerator.* sinv;

% Fix up dp == 1 which causes NaN quaternions. This means the two
% quaternions are the same - just use the first.
%
% If sinv is inf, qo is nan. But we can't look for nans in qo because they
% may have come in from q1 or q2 == nan at input. Instead find infs in
% sinv, then expand it to match size of qo, then put the expanded q1s where
% we found infs.
infmap = isinf(sinv);
if any(infmap(:)) % Don't do this unless necessary
    infmapExpanded = infmap & true(size(qo,1)); % handle implicit expansion to size(qo)
    infmapExpandedNumeric = cast(infmapExpanded, 'double'); % same thing as above with 1s
    replaceval = q1 .* infmapExpandedNumeric; % an array with q1 everywhere there's an inf in infmapExpanded
    qo(infmapExpanded,:) = replaceval(infmapExpanded,:); % replace
end

qo = my_quaternion.quat_normalize(qo);

end