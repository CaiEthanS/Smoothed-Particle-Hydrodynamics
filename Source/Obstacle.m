function [] = Obstacle(xMax, yMax, h, N, kappa, mu, rho_0, beta, dt)
% 2D Smoothed Particle Hydrodynamics simulation for a block of fluid of
% stiffness kappa, viscosity mu, rest density rho_0, and boundary damping
% beta, positioned to be above a circle centered at [0.4, 1]. The particles
% are contained in a grid of size xMax by yMax.

% Create video object
% vid = VideoWriter('ObjectCase', 'MPEG-4');
% vid.FrameRate = 30;
% vid.Quality = 100;
% open(vid);

% Other initial condition (besides the ones inputted)
m = rho_0 / N;
boundParticles = 50;

% Set grid size
Nx = floor(xMax / h);
Ny = floor(yMax / h);
dx = xMax / Nx;
dy = yMax / Ny;

particles(1:N + boundParticles) = struct('pos', [], 'vel', [0, 0], 'force', [0, 0], 'rho', [], 'neighbors', [], 'neighborsDist', []);
sp = 0.09;
k = 1;
while k <= N
    for r = 1:ceil(sqrt(N))
        for c = 1:ceil(sqrt(N))
            particles(k).pos = [r * sp + 0.5, c * sp + 0.7];
            k = k + 1;
        end
    end
end
theta = 0;
for f = 1:boundParticles
    particles(f + N).pos = [0.25*cos(theta) + 1, 0.25*sin(theta) + 0.4];
    theta = theta + 2*pi/boundParticles;
end
N = N + boundParticles;
% Calculate the bin arrangement (.adjacentBins)
numBins = Nx * Ny;
bins(1:numBins) = struct('particleIDs', [], 'adjacentBins', []);
for z = 1:numBins
    bins(z).adjacentBins = getAdjacentBins(z, numBins, Nx, Ny);
end

% for all timesteps
for t = 0:dt:3
    % Hash particles to build the bins.particleIDs lists
    for k = 1:N
        binNum = (ceil(particles(k).pos(1)/dx) - 1) * Ny + ceil((yMax - particles(k).pos(2))/dy);
        bins(binNum).particleIDs = [bins(binNum).particleIDs, k];
    end
    % Build the list of neighbors for each particle
    for z = 1:numBins
        % Get all particles in bin z;
        AList = bins(z).particleIDs;
        % Get vector of bin IDs adjacent to bin z;
        adjacentBins = bins(z).adjacentBins;
        w = [z, adjacentBins];
        consider = [];
        for a = 1:length(w)
            % Get all particles in bin w;
            consider = [consider, bins(w(a)).particleIDs];
        end
        for k = 1:length(AList) % for each particle in bin z
            for j = 1:length(consider) % for each possible neighboring particle
                xDiff = particles(AList(k)).pos(1) - particles(consider(j)).pos(1);
                yDiff = particles(AList(k)).pos(2) - particles(consider(j)).pos(2);
                dist = sqrt(xDiff^2 + yDiff^2);
                if dist < h && AList(k) ~= consider(j)
                    % Particle j is a neighbor of k, store distance for rho
                    particles(AList(k)).neighbors = [particles(AList(k)).neighbors, consider(j)];
                    particles(AList(k)).neighborsDist = [particles(AList(k)).neighborsDist ,dist];
                end
            end
        end
    end
    % Density calculation
    for k = 1:N
        particles(k).rho = 4*m/(pi*h^2);
        for n = 1:length(particles(k).neighbors)
            % Summation of neighboring particle density contributions
            neighborRho = (h^2 - particles(k).neighborsDist(n)^2)^3;
            particles(k).rho = particles(k).rho + (4*m/(pi*h^8))*neighborRho;
        end
    end
    % Force calculation
    for k = 1:N
        rhoK = particles(k).rho;
        for n = 1:length(particles(k).neighbors)
            %  Calculate the total force on each particle
            j = particles(k).neighbors(n);
            q = particles(k).neighborsDist(n)/h;
            rhoJ = particles(j).rho;
            posDiff = particles(k).pos - particles(j).pos;
            velDiff = particles(k).vel - particles(j).vel;
            Fkj = m*(1 - q)*(15*kappa*(rhoK + rhoJ - 2*rho_0)*(1 - q)*posDiff/q - 40*mu*velDiff)/(pi*(h^4)*rhoJ);
            particles(k).force = particles(k).force + Fkj;
        end
        FExt = [0, -9.8]*rhoK;
        particles(k).force = particles(k).force + FExt;
    end
    for k = 1:N - boundParticles
        %  Update the velocity and position of each particle
        particles(k).vel = (particles(k).force/particles(k).rho)*dt + particles(k).vel;
        particles(k).pos = particles(k).vel * dt + particles(k).pos;
        %  Apply any necessary boundary conditions
        xPos = particles(k).pos(1);
        yPos = particles(k).pos(2);
        if xPos < 0
            particles(k).pos(1) = -xPos;
            particles(k).vel(1) = -beta * particles(k).vel(1);
        elseif xPos > xMax
            particles(k).pos(1) = 2*xMax - xPos;
            particles(k).vel(1) = -beta * particles(k).vel(1);
        end
        if yPos < 0
            particles(k).pos(2) = -yPos;
            particles(k).vel(2) = -beta * particles(k).vel(2);
        elseif yPos > yMax
            particles(k).pos(2) = 2*yMax - yPos;
            particles(k).vel(2) = -beta * particles(k).vel(2);
        end
    end
    %  Visualize the particle system and save frame into video
    x = zeros(1, N);
    y = zeros(1, N);
    for k = 1:N
        x(k) = particles(k).pos(1);
        y(k) = particles(k).pos(2);
    end
    drawnow
    plot(x(1:N - boundParticles), y(1:N - boundParticles), 'bo', x(N - (boundParticles - 1): N), y(N - (boundParticles - 1): N), 'ko');
    grid on
    xlim([0,xMax]);
    ylim([0,yMax]);
%     writeVideo(vid, getframe(gcf));
    %  Empty bins and neighbors for next time iteration
    for z = 1:numBins
        bins(z).particleIDs = [];
    end
    for k = 1:N
        particles(k).force = [0, 0];
        particles(k).neighbors = [];
        particles(k).neighborsDist = [];
    end
end
% close(vid);
end