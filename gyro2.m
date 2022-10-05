% Gyroscope Simulation

clear all
close all

g = 9.8;														% Gravitational acceleration (m/s^2)
omega = 10;														% Initial angular velocity (rad/s)

r = 1;															% Radius of wheel (m)
a = 0.5;														% Length of axle (m)
n = 20;															% Number of points on the rim
M_rim = 1;														% Total mass of the rim (kg)
M_axle = 1;														% Total mass of the axle (kg)
S_rim = 2000 * (M_rim / n) * omega ^ 2;							% Stiffness of each rim link (kg/s^2)
D_rim = 0 * (M_rim / n) * omega;								% Damping constant of each rim link (kg/s)
S_spoke = S_rim;												% Stiffnes of each spoke (kg/s^2)
D_spoke = D_rim;												% Damping constant of each spoke (kg/s)
S_axle =  S_rim;												% Stiffness of the axle (kg/s^2)
D_axle =  D_rim;												% Damping constant of the axle (kg/s)

[kmax, lmax, X, jj, kk, S, D, Rzero, M] = ...
	wheel(r, a, n, M_rim, M_axle, S_rim, D_rim, S_spoke, D_spoke, S_axle, D_axle);

figure(1)														% Setup for animation
lr = 1 : n;														% Rim links
ls = n + 1 : 3 * n;												% Spoke links
la = 3 * n + 1;													% Axle link
nskip = 5;

hr = plot3([X(jj(lr), 1), X(kk(lr), 1)]', [X(jj(lr), 2), X(kk(lr), 2)]', ...
	[X(jj(lr), 3), X(kk(lr), 3)]', "r", "linewidth", 3);
axis([-1.2 * r, 1.2 * r, -1.2 * r, 1.2 * r, -1.2 * a, 1.2 * a]);
hold on
hs = plot3([X(jj(ls), 1), X(kk(ls), 1)]', [X(jj(ls), 2), X(kk(ls), 2)]', ...
	[X(jj(ls), 3), X(kk(ls), 3)]', "b");
ha = plot3([X(jj(la), 1), X(kk(la), 1)]', [X(jj(la), 2), X(kk(la), 2)]', ...
	[X(jj(la), 3), X(kk(la), 3)]', "k", "linewidth", 3);
hold off
axis equal
axis manual
axis(1.2 * [-r, r, -r, r, -a, a])
drawnow


fixed = n + 1;													% Index of point that will be held fixed
forced = n + 2;													% Index of point to which external force will be applied

tmax = 20;														% Duration of simulation (s)
clockmax = 10000;												% Number of time steps
dt = tmax / clockmax;											% (s)
X_save = zeros(clockmax, 3);
t_save = zeros(clockmax, 1);

t_ext_start = 2;												% Time that external force starts (s)
t_ext_stop =  3;												% Time that external force stops (s)
F_ext = 0.5 * [1, 0, 0] * M_rim * g;							% External force (kg·m/s^2)
U = omega * [-X(:, 2), X(:, 1), zeros(kmax, 1)];				% Initial velocity

for clock = 1 : clockmax

	t = clock * dt;
	DX = X(jj, :) - X(kk, :);									% Link vectors
	DU = U(jj, :) - U(kk, :);									% Link velocity difference vectors
	R = sqrt(sum(DX .^ 2, 2));									% Link lengths
  	T = S .* (R - Rzero) + (D ./ R) .* sum(DX .* DU, 2);		% Link tensions
	TR = T ./ R;												% Link tensions divided by link lengths
	FF = [TR, TR, TR] .* DX;									% Link force vectors

	F = zeros(kmax, 3);											% Initialize force array for mass points 
	F(:, 3) = -M * g;											% Apply gravitational force to each link point

	for link = 1 : lmax
		F(kk(link), :) = F(kk(link), :) + FF(link, :);			% Add force of the link to one end
		F(jj(link), :) = F(jj(link), :) - FF(link, :);			% Add force of the link to the other end
	end

	if (t > t_ext_start) && (t < t_ext_stop)
		F(forced, :) = F(forced, :) + F_ext;					% Apply external force during specified time interval
	end

	U = U + dt * F ./ [M, M, M];								% Update velocities of all points
	U(fixed, :) = 0; 											% The fixed point must have velocity zero
	X = X + dt * U;												% Update positions of all points
	X_save(clock, :) = X(n + 2, :);								% Store positions for future plotting
	t_save(clock) = t;											% Store current time for future plotting
  
	%%%%%%%%%% %%%%%%%%%% %%%%%%%%%% %%%%%%%%%% %%%%%%%%%% %%%%%%%%%% %%%%%%%%%% %%%%%%%%%% %%%%%%%%%% %%%%%%%%%%
	%%                                                Animation                                                %%
	%%%%%%%%%% %%%%%%%%%% %%%%%%%%%% %%%%%%%%%% %%%%%%%%%% %%%%%%%%%% %%%%%%%%%% %%%%%%%%%% %%%%%%%%%% %%%%%%%%%%
	if mod(clock, nskip) == 0
		c = 0;
		for l = lr
			c = c + 1;
			hr(c).XData = [X(jj(l), 1), X(kk(l), 1)];
			hr(c).YData = [X(jj(l), 2), X(kk(l), 2)];
			hr(c).ZData = [X(jj(l), 3), X(kk(l), 3)];
		end
		c = 0;
		for l = ls
			c = c + 1;
			hs(c).XData = [X(jj(l), 1), X(kk(l), 1)];
			hs(c).YData = [X(jj(l), 2), X(kk(l), 2)];
			hs(c).ZData = [X(jj(l), 3), X(kk(l), 3)];
		end
		c = 0;
		for l = la
			c = c + 1;
			ha(c).XData = [X(jj(l), 1), X(kk(l), 1)];
			ha(c).YData = [X(jj(l), 2), X(kk(l), 2)];
			ha(c).ZData = [X(jj(l), 3), X(kk(l), 3)];
		end
		drawnow
    end

end

figure(2)
plot(t_save', X_save')