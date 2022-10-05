% This function constructs a wheel with n points on the rim and 2 axle points, one on each side of the plane of the rim.
% Each rim point is linked to its two neighbors on the rim, and it is also linked to each axle point by a spoke.
% Finally, the two axle points are linked to each other.

function [kmax, lmax, X, jj, kk, S, D, Rzero, M] = ...
	wheel(r, a, n, M_rim, M_axle, S_rim, D_rim, S_spoke, D_spoke, S_axle, D_axle)

% Input parameters:
%	r															Radius of wheel
%	a															Length of axle
%	n															Number of points on the rim
%	M_rim														Total mass of the rim
%	M_axle														Total mass of the axle
%	S_rim														Stiffness of each rim link
%	D_rim														Damping constant of each rim link
%	S_spoke														Stiffnes of each spoke
%	D_spoke														Damping constant of each spoke
%	S_axle														Stiffness of the axle
%	D_axle														Damping constant of the axle

% Output values:
%	kmax														Number of points
%	lmax														Number of links
%	X(k, :)														Coordinates of point k
%	jj(l), kk(l)												Indices of points connected by link l
%	S(l)														Stiffness of link l
%	D(l)														Damping constant of link l
%	Rzero(l)													Rest length of link l
%	M(k)														Mass of point k

kmax = n + 2;													% Rim points n and axle points 2
lmax = 3 * n + 1;												% Rim links n, spokes 2n, and axle link 1
R_spoke = sqrt(r ^ 2 + (a / 2) ^ 2);							% Length of each spoke
R_rim = 2 * r * sin((2 * pi) / (2 * n))							% Length of each rim link (chord length)

X = zeros(n + 2, 3);											% Initialization of coordinates
for k = 1 : n
	theta = 2 * pi * k / n;
	X(k, :) = r * [cos(theta), sin(theta), 0];					% Coordinates of points on the rim
end
X(n + 1, :) = (a / 2) * [0, 0, -1];								% Coordinate of one end of the axle
X(n + 2, :) = (a / 2) * [0, 0, 1];								% Coordinate of the other end of the axle

jj = zeros(lmax, 1); kk = zeros(lmax, 1);						% Initialization of links
for k = 1 : n
	jj(k) = k; kk(k) = k + 1 - (k == n) * n;					% Link to next point on rim    
	jj(n + k) = k; kk(n + k) = n + 1;							% Link to one axle point
	jj(2 * n + k) = k; kk(2 * n + k) = n + 2;					% Link to other axle point
end
jj(3 * n + 1) = n + 1; kk(3 * n + 1) = n + 2;					% Link between axle points

M = [(M_rim / n) * ones(n, 1); (M_axle / 2) * ones(2, 1)];		% Mass of each point
S = [S_rim * ones(n, 1); S_spoke * ones(2 * n, 1); S_axle];		% Stiffness of each link
D = [D_rim * ones(n, 1); D_spoke * ones(2 * n, 1); D_axle];		% Damping constant of each link
Rzero = [R_rim * ones(n, 1); R_spoke * ones(2 * n, 1); a];		% Rest length of each link
