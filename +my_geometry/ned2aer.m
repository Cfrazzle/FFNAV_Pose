function [az, elev, slantRange] = ned2aer(xNorth, yEast, zDown, angleUnit)
%NED2AER Local Cartesian NED to spherical AER
%
%   [az, elev, slantRange] = NED2AER(xNorth, yEast, zDown) transforms point
%   locations in 3-D from local Cartesian coordinates (xNorth, yEast,
%   zDown) to local spherical coordinates (azimuth angle, elevation angle,
%   and slant range). The output angles are returned in degrees.
%
%   [...] = NED2AER(..., angleUnit) uses angleUnit, which matches either
%   'degrees' or 'radians', to specify the units of the azimuth and
%   elevation angles.
%
%   North-east-down (NED) is a right-handed local Cartesian system with the
%   X-axis directed north and parallel to the local tangent plane, the
%   Y-axis directed east, and the Z-axis directed downward along the local
%   normal to the ellipsoid.
%
%   As always, the azimuth angle is measured clockwise (east) from north,
%   from the perspective of a viewer looking down on the local horizontal
%   plane. Equivalently, in the case of NED, it is measured clockwise from
%   the positive X-axis in the direction of the positive Y-axis. The
%   elevation angle is the angle between a vector from the origin and the
%   local horizontal plane. The slant range is the 3-D Euclidean distance
%   from the origin.
%
%   The transformation is similar to CART2SPH, except that the output
%   angles are in degrees by default, the Z-axis is directed downward, and
%   the first output is in the half-open interval [0 360) in degrees, or
%   [0 2*pi) in radians.
%
%   Class support for inputs xNorth, yEast, zDown:
%      float: double, single
%
%   See also AER2NED, ENU2AER, NED2ECEF, NED2GEODETIC

% Copyright 2012-2020 The MathWorks, Inc.

%#codegen

if nargin < 4 || map.geodesy.isDegree(angleUnit)
    atan2fun = @atan2d;
else
    atan2fun = @atan2;
end

[az, elev, slantRange] = enu2aerFormula(yEast, xNorth, -zDown, atan2fun);

end