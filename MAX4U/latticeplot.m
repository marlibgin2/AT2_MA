function latticeplot(RING) 
% It plots the lattice function along with magnet blocks & their name
%% Inputs 
% RING: Lattice structure of an achromat
%% Usage examples
% latticeplot(RING)

%% History
% Saroj 2024/11/12 

% Define magnet families and their properties
height = 0.6;
magnetFamilies = struct( ...
    'Q1', struct('Color', 'r',       'Height', height), ...
    'Q2', struct('Color', 'r',       'Height', height), ...
    'Q3', struct('Color', 'r',       'Height', height), ...
    'Q4', struct('Color', 'r',       'Height', height), ...
    'D1', struct('Color', '#0072BD', 'Height', height), ...
    'D2', struct('Color', '#0072BD', 'Height', height), ...
    'D3', struct('Color', '#0072BD', 'Height', height), ...
    'R1', struct('Color', '#EDB120', 'Height', height), ...
    'R2', struct('Color', '#EDB120', 'Height', height), ...
    'R3', struct('Color', '#EDB120', 'Height', height), ...
    'S1', struct('Color', '#77AC30', 'Height', height), ...
    'S2', struct('Color', '#77AC30', 'Height', height), ...
    'S3', struct('Color', '#77AC30', 'Height', height), ...
    'S4', struct('Color', '#77AC30', 'Height', height), ...
    'S5', struct('Color', '#77AC30', 'Height', height), ...
    'S6', struct('Color', '#77AC30', 'Height', height), ...
    'O1', struct('Color', 'm',       'Height', height), ...
    'O2', struct('Color', 'm',       'Height', height), ...
    'O3', struct('Color', 'm',       'Height', height));

% Iterate over each family in the structure
familyNames = fieldnames(magnetFamilies);
hold on;

for f = 1:numel(familyNames)
    familyName = familyNames{f};
    familyProps = magnetFamilies.(familyName);

    % Find elements of this family in RING
    familyIndices = findcells(RING, 'FamName', familyName);

    % Group contiguous elements for labeling
    if ~isempty(familyIndices)
        groups = splitIntoContiguousGroups(familyIndices);
        for g = 1:numel(groups)
            groupIndices = groups{g};
            startIdx = groupIndices(1);
            endIdx = groupIndices(end);
            
            % Calculate the position and length of this contiguous block
            startPos = findspos(RING, startIdx);
            endPos = findspos(RING, endIdx) + RING{endIdx}.Length;
            blockLength = endPos - startPos;
            blockCenter = startPos + blockLength / 2;

            % Plot each magnet in the contiguous block
            for i = groupIndices
                elemPos = findspos(RING, i);
                elemLength = RING{i}.Length;
                rectangle('Position', [elemPos 0 elemLength familyProps.Height], ...
                    'EdgeColor', familyProps.Color, 'FaceColor', familyProps.Color);
            end

            % Add text label at the center of the block
            text(blockCenter, familyProps.Height / 2, familyName, ...
                'fontsize', 10, 'fontweight', 'bold', 'Rotation', 90, ...
                'HorizontalAlignment', 'center');
        end
    end
end

% Calculate and plot the Twiss parameters
RING=splitlat(RING,10);
[TD, ~] = twissring(RING, 0, 1:length(RING)+1, 'chrom');
sPos = cat(1, TD.SPos);
beta = cat(1, TD.beta);
dispersion = cat(2, TD.Dispersion);

% Plot beta functions
yyaxis left
plot(sPos, beta(:,1), 'b-', sPos, beta(:,2), 'r-', 'LineWidth', 2);
ylabel('\beta [m]');

% Plot horizontal dispersion
yyaxis right
plot(sPos, dispersion(1,:), 'k:', 'LineWidth', 2);
xlabel('Path length [m]');
ylabel('\eta [m]');
legend('\beta_x', '\beta_y', '\eta_x');
set(gca, 'FontSize', 14);
xlim([0, max(sPos)]);

hold off;
end

function groups = splitIntoContiguousGroups(indices)
% This function takes a list of indices and splits them into contiguous groups.
    groups = {};
    currentGroup = indices(1);
    
    for i = 2:length(indices)
        if indices(i) == indices(i-1) + 1
            currentGroup = [currentGroup, indices(i)];
        else
            groups{end+1} = currentGroup; 
            currentGroup = indices(i);
        end
    end
    groups{end+1} = currentGroup;
end
