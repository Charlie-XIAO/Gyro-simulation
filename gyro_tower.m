%%%%%%%%%% %%%%%%%%%% %%%%%%%%%% %%%%%%%%%% %%%%%%%%%% %%%%%%%%%% %%%%%%%%%% %%%%%%%%%% %%%%%%%%%% %%%%%%%%%%
%%                                       Gyroscope Tower Simulation                                        %%
%%%%%%%%%% %%%%%%%%%% %%%%%%%%%% %%%%%%%%%% %%%%%%%%%% %%%%%%%%%% %%%%%%%%%% %%%%%%%%%% %%%%%%%%%% %%%%%%%%%%

clearvars
clear global
close all

ngyro = input("How many gyroscopes in the gyro-tower? [int >= 1] ");
while mod(ngyro, 1) ~= 0 || ngyro < 1
    ngyro = input("Invalid. How many gyroscopes in the gyro-tower? [int >= 1] ");
end
reverse = input("Do you want each gyroscope to spin in reverse direction? [true/false] ");
while reverse ~= true && reverse ~= false
    reverse = input("Invalid. Do you want to save a movie? [true/false] ");
end
video = input("Do you want to save a movie? [true/false] ");
while video ~= true && video ~= false
    video = input("Invalid. Do you want to save a movie? [true/false] ");
end
if video
    if reverse
        writerObj = VideoWriter(sprintf("./videos/gyro-tower-%d-rev.mp4", ngyro), "MPEG-4");
    else
        writerObj = VideoWriter(sprintf("./videos/gyro-tower-%d-irr.mp4", ngyro), "MPEG-4");
    end
    writerObj.FrameRate = 30;
    open(writerObj);
end

g = 9.8;														% Gravitational acceleration (m/s^2)
omega = 10;                                                     % Initial angular velocity (rad/s)
r = 1;                                                          % Radius of wheel (m)
a = 0.5;                                                        % Length of axle (m)
n = 20;                                                         % Number of points on the rim
M_rim = 1;                                                      % Total mass of the rim (kg)
M_axle = 1;                                                     % Total mass of the axle (kg)
S_rim = 2000 * (M_rim / n) * omega ^ 2;                         % Stiffness of each rim link (kg/s^2)
D_rim = 0 * (M_rim / n) * omega;                                % Damping constant of each rim link (kg/s)
S_spoke = S_rim;                                                % Stiffnes of each spoke (kg/s^2)
D_spoke = D_rim;                                                % Damping constant of each spoke (kg/s)
S_axle =  S_rim;                                                % Stiffness of the axle (kg/s^2)
D_axle =  D_rim;                                                % Damping constant of the axle (kg/s)
SS = S_rim * 10;                                                % Stiffness of the invisible link for simplicity (kg/s^2)
DD = D_rim * 10;                                                % Damping constant of the invisible link for simplicity (kg/s)
RRzero = a / 1000;                                              % Rest length of the invisible link for simplicity (m)

figure(1)                                                       % Setup for animation
nskip = 5;                                                      % Clock skip
lr = zeros(1, ngyro * n);                                       % Rim links
ls = zeros(1, ngyro * 2 * n);                                   % Spoke links
la = zeros(1, ngyro);                                           % Axle link
kmax = 0;                                                       % Number of points in total
lmax = 0;                                                       % Number of links in total
X = zeros((n + 2) * ngyro, 3);                                  % Coordinates of points
jj = zeros((3 * n + 1) * ngyro, 1);                             % Indices of points on one end of the links
kk = zeros((3 * n + 1) * ngyro, 1);                             % Indices of points on the other end of the links
S = zeros((3 * n + 1) * ngyro, 1);                              % Stiffness of links
D = zeros((3 * n + 1) * ngyro, 1);                              % Damping constant of links
Rzero = zeros((3 * n + 1) * ngyro, 1);                          % Rest length of links
M = zeros((n + 2) * ngyro, 1);                                  % Mass of points
for gyro = 1 : ngyro                                            % Update information for each gyroscope
    lr((gyro - 1) * n + 1 : gyro * n) = (gyro - 1) * (3 * n + 1) + 1 : (gyro - 1) * (3 * n + 1) + n;
    ls((gyro - 1) * 2 * n + 1 : gyro * 2 * n) = (gyro - 1) * (3 * n + 1) + n + 1 : (gyro - 1) * (3 * n + 1) + 3 * n;
    la(gyro) = (gyro - 1) * (3 * n + 1) + 3 * n + 1;
    [ksing, lsing, ...                                                  % Number of points and links within a single gyroscope 
        X((gyro - 1) * (n + 2) + 1 : gyro * (n + 2), :), ...            % Coordinates of points within a single gyroscope
        jj((gyro - 1) * (3 * n + 1) + 1 : gyro * (3 * n + 1)), ...      % Indices of points on one end of the links within a single gyroscope
        kk((gyro - 1) * (3 * n + 1) + 1 : gyro * (3 * n + 1)), ...      % Indices of points on the other end of the links within a single gyroscope
        S((gyro - 1) * (3 * n + 1) + 1 : gyro * (3 * n + 1)), ...       % Stiffness of links within a single gyroscope
        D((gyro - 1) * (3 * n + 1) + 1 : gyro * (3 * n + 1)), ...       % Damping constant of links within a single gyroscope
        Rzero((gyro - 1) * (3 * n + 1) + 1 : gyro * (3 * n + 1)), ...   % Rest length of links within a single gyroscope
        M((gyro - 1) * (n + 2) + 1 : gyro * (n + 2))] = ...             % Mass of points within a single gyroscope
        wheel(r, a, (gyro - 1) * (a + RRzero), n, M_rim, M_axle, S_rim, D_rim, S_spoke, D_spoke, S_axle, D_axle, (gyro - 1) * (n + 2));
    kmax = kmax + ksing;                                                % Update total number of points
    lmax = lmax + lsing;                                                % Update total number of links
end

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
axis(1.2 * [-r, r, -r, r, -0.5, a * ngyro])
drawnow

tmax = 10;                                                      % Duration of simulation (s)
clockmax = 50000;                                               % Number of time steps
dt = tmax / clockmax;                                           % Time step (s)
t_save = zeros(clockmax, 1);                                    % Save each time step
KE_save = zeros(clockmax, 1);                                   % Save kinetic energy for each time step
GPE_save = zeros(clockmax, 1);                                  % Save gravitational potential energy for each time step
EPE_save = zeros(clockmax, 1);                                  % Save elastic potential energy for each time step
E_save = zeros(clockmax, 1);                                    % Save total energy for each time step

t_ext_start = 1;                                                % Time that external force starts (s)
t_ext_stop = 2;                                                 % Time that external force stops (s)
F_ext = 0.5 * [1, 0, 0] * M_rim * g;                            % External force (kg·m/s^2)
forced = ngyro * (n + 2);                                       % The point where external force is applied
if reverse
    U = zeros(kmax, 3);                                                                 % Initial velocity
    for gyro = 1 : ngyro
        tmptmp = (gyro - 1) * ksing + 1 : gyro * ksing;                                 % Temporary list for update
        if mod(gyro, 2) == 0
            U(tmptmp, :) = omega * [-X(tmptmp, 2), X(tmptmp, 1), zeros(ksing, 1)];      % Update in one direction
        else
            U(tmptmp, :) = omega * [X(tmptmp, 2), -X(tmptmp, 1), zeros(ksing, 1)];      % Update in the reverse direction
        end
    end
else
    U = omega * [-X(:, 2), X(:, 1), zeros(kmax, 1)];                                    % Initial velocity all same direction
end

for clock = 1 : clockmax

    t = clock * dt;                                             % Current time
    DX = X(jj, :) - X(kk, :);                                   % Link vectors
    DU = U(jj, :) - U(kk, :);                                   % Link velocity difference vectors
    R = sqrt(sum(DX .^ 2, 2));                                  % Link lengths
    DR = R - Rzero;                                             % Elastic changes of lengths
    T = S .* DR + (D ./ R) .* sum(DX .* DU, 2);                 % Link tensions
    TR = T ./ R;                                                % Link tensions divided by link lengths
    FF = [TR, TR, TR] .* DX;                                    % Link force vectors

    F = zeros(kmax, 3);                                         % Initialize force array for mass points 
    F(:, 3) = -M * g;                                           % Apply gravitational force to each link point

    for link = 1 : lmax
        F(kk(link), :) = F(kk(link), :) + FF(link, :);          % Add force of the link to one end
        F(jj(link), :) = F(jj(link), :) - FF(link, :);          % Add force of the link to the other end
    end

    if (t > t_ext_start) && (t < t_ext_stop)                    % In the specified time interval
        F(forced, :) = F(forced, :) + F_ext;                    % Apply external force during specified time interval
    end

    EE = 0;                                                                 % Total elastic potential energy within invisible links
    for gyro = 2 : ngyro                                                    % Dealing with invisible links between top of the lower gyroscope and bottom of the upper one
        DXX = X(gyro * (n + 2) - 1, :) - X((gyro - 1) * (n + 2), :);        % Link vector
        DUU = U(gyro * (n + 2) - 1, :) - U((gyro - 1) * (n + 2), :);        % Link velocity difference vector
        RR = sqrt(sum(DXX .^ 2, 2));                                        % Link length
        DRR = RR - RRzero;                                                  % Elastic change of length
        TT = SS * DRR + (DD / RR) * sum(DXX .* DUU, 2);                     % Link tension
        TTRR = TT / RR;                                                     % Link tension divided by link length
        FFFF = [TTRR, TTRR, TTRR] .* DXX;                                   % Link force vector
        F(gyro * (n + 2) - 1, :) = F(gyro * (n + 2) - 1, :) - FFFF;         % Add force of the link to one end
        F((gyro - 1) * (n + 2), :) = F((gyro - 1) * (n + 2), :) + FFFF;     % Add force of the link to the other end
        EE = EE + 1 / 2 * SS * DRR ^ 2;                                     % Add elastic potential energy to the total potential energy within invisible links
    end

    U = U + dt * F ./ [M, M, M];                                % Update velocities of all points
    U(n + 1, :) = 0;                                            % The fixed point must have velocity zero
    X = X + dt * U;                                             % Update positions of all points

    t_save(clock) = t;                                          % Save current time step
    E_save(clock) = ...                                         % Save current total energy
        1 / 2 * sum(M .* sqrt(sum(U .^ 2, 2)) .^ 2) + ...       % Total kinetic energy
        sum(M .* X(:, 3) * g) + ...                             % Total gravitational potential energy
        1 / 2 * sum(S .* DR .^ 2) + EE;                         % Total elastic potential energy
  
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

        if video
            frame = getframe(gcf);
            writeVideo(writerObj, frame);
        end

    end

end

figure(2)
plot(t_save', E_save')

if video
    close(writerObj);
end
