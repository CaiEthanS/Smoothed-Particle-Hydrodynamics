function [adjacentBins] = getAdjacentBins(binNum, binNumMax, M, N)
% Returns the adjacentBins for a certain bin in a spatially hashed region
% with MxN number of bins total.
adjacentBins = zeros(1, 8);

% Left wall
if binNum <= M
    adjacentBins(1) = binNum + M;
% Top left corner
    if binNum == 1
        adjacentBins(2) = binNum + 1;
        adjacentBins(3) = binNum + M + 1;
    % Bottom left corner    
    elseif binNum == M
        adjacentBins(2) = binNum - 1;
        adjacentBins(3) = binNum + M - 1;
    % Middle of left wall    
    else    
        adjacentBins(2) = binNum - 1;
        adjacentBins(3) = binNum + 1;
        adjacentBins(4) = binNum + M - 1;
        adjacentBins(5) = binNum + M + 1;
    end
% Right wall
elseif binNum > binNumMax - M
    adjacentBins(1) = binNum - M;
    if binNum == binNumMax - M + 1
        % Right top corner
        adjacentBins(2) = binNum + 1;
        adjacentBins(3) = binNum - M + 1;
    elseif binNum == binNumMax
        % Right bottom corner
        adjacentBins(2) = binNum - 1;
        adjacentBins(3) = binNum - M - 1;
    else
        % Middle of right wall
        adjacentBins(2) = binNum - 1;
        adjacentBins(3) = binNum + 1;
        adjacentBins(4) = binNum - M - 1;
        adjacentBins(5) = binNum - M + 1;
    end
% Top Wall
elseif mod(binNum, M) == 1
    adjacentBins(1) = binNum + 1;
    % Top left corner
    if binNum == 1
        adjacentBins(2) = binNum + M;
        adjacentBins(3) = binNum + M + 1;
    % Top right corner
    elseif binNum == (N - 1) * M + 1
        adjacentBins(2) = binNum - M;
        adjacentBins(3) = binNum - M + 1;
    % Middle of top wall
    else
        adjacentBins(2) = binNum + M + 1;
        adjacentBins(3) = binNum + M - 2;
        adjacentBins(4) = binNum + M;
        adjacentBins(5) = binNum - M;
    end
% Bottom Wall
elseif mod(binNum, M) == 0
    adjacentBins(1) = binNum - 1;
    % Bottom left corner
    if binNum == M
        adjacentBins(2) = binNum + M - 1;
        adjacentBins(3) = binNum + M;
    % Bottom right corner
    elseif binNum == binNumMax
        adjacentBins(2) = binNum - M;
        adjacentBins(3) = binNum - M - 1;
    % Middle of bottom wall
    else
        adjacentBins(2) = binNum - M - 1;
        adjacentBins(3) = binNum - M;
        adjacentBins(4) = binNum + M - 1;
        adjacentBins(5) = binNum + M;
    end
% Middle of grid, not on any walls
else
    adjacentBins(1) = binNum - M - 1;
    adjacentBins(2) = binNum - M;
    adjacentBins(3) = binNum - M + 1;
    adjacentBins(4) = binNum - 1;
    adjacentBins(5) = binNum + 1;
    adjacentBins(6) = binNum + M - 1;
    adjacentBins(7) = binNum + M;
    adjacentBins(8) = binNum + M + 1;
end

adjacentBins = sort(nonzeros(adjacentBins)).';
end