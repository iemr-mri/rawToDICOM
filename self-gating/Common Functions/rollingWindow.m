function rollingWindow(rawCardiac, filteredCardiac, rawBreath, filteredBreath, roughPeakList, finePeakList, breathStarts, timeline, pm, showRaw)
% Takes in curve before and after filtering and and rolls through it 
% section by section.

disp('Section 6: rollingWindow');

fullAxis            = 1:length(rawCardiac);
secsPerWindow       = 2.5;                                                  % How many seconds should each window show
windowPerStep       = 0.02;                                                 % How much of the window we move per step
waitTime            = 0.1;                                                  % How long to wait between steps
timepointsPerWindow = round(secsPerWindow * 1000 / pm.newTempRes);          % Convert to milliseconds 
stepSize            = round(timepointsPerWindow * windowPerStep);
numSteps            = floor(length(filteredCardiac) / stepSize) - 1;        % Take off last step to not overflow the graph
graphScaler         = max(rawCardiac) / max(filteredCardiac);               % Make graphs equalish height. Quite violent, but should work for stable graphs
filteredCardiac     = filteredCardiac * graphScaler;
cardiacBreathScale  = (max(filteredCardiac) - min(filteredCardiac)) / (max(filteredBreath) - min(filteredBreath));
filteredBreath      = filteredBreath * cardiacBreathScale;
rawBreath           = rawBreath * cardiacBreathScale;
breathStarts(2,:)   = breathStarts(2,:) * cardiacBreathScale;
lowest              = min([min(filteredCardiac) min(filteredBreath)]);
highest             = max([max(filteredCardiac) max(filteredBreath)]);
ylimits             = [lowest,highest]*1.05;  % Shady but works well for stable curve that oscillates around the same values

gifFileName = [pm.animationFolder pm.expName '/' pm.scanNumber ' rolling.gif'];

% Remove pre-existing GIF
if exist(gifFileName, 'file') == 2
    delete(gifFileName);
    disp([gifFileName ' deleted for new GIF.'])
end

% Setup figure
figHandle = figure('Position', [100,100,1800,800]);
title('PCA of cardiac motion');
xlabel('Time [ms]');
ylabel('PC Value [a.u.]');
grid on;
hold on;

% Create stop button and set its callback % K: Pure GPT, don't know exactly
% how this GUI stuff works!
hButton = uicontrol('Style', 'pushbutton', 'String', 'Stop Animation', ...
                    'Position', [20, 350, 150, 125], ...
                    'Callback', @stopAnimation);
stopFlag = false;

% Nested function to handle the button press
function stopAnimation(~, ~)
    stopFlag = true;
end

% Loop through the graph in chunks
for step = 1:numSteps-floor(timepointsPerWindow/stepSize)

    if stopFlag == true
        disp(['Animation stopped.' newline]);
        return;
    end

    timeStart               = step*stepSize;
    timeEnd                 = timeStart+timepointsPerWindow;
    timeRange               = timeStart:timeEnd;
    timeChunk               = timeline(timeRange);
    rawCardiacChunk         = rawCardiac(timeRange);
    filteredCardiacChunk    = filteredCardiac(timeRange);
    rawBreathChunk          = rawBreath(timeRange);
    filteredBreathChunk     = filteredBreath(timeRange);
    
    % Clear axes to make room for new timestep
    cla;

    % Plot current chunk
    hold on;
    plot(timeChunk,filteredCardiacChunk, '-', 'LineWidth',6, 'Color','#ff0065');
    plot(timeChunk,filteredBreathChunk, '-', 'LineWidth',5, 'Color','#0096FF');

    if showRaw == true
        plot(timeChunk,rawCardiacChunk, '-', 'LineWidth',2, 'Color','#a72707');
        plot(timeChunk,rawBreathChunk, '-', 'LineWidth',2, 'Color','#ADD8E6');
    end

    % % Detect and mark rough peaks in the current chunk
    % roughInRange    = (roughPeakList(:,1) >= timeStart) & (roughPeakList(:,1) <= timeEnd);
    % 
    % if any(roughInRange)
    %     roughPeakPos = roughPeakList(roughInRange,:);
    %     for peak = 1:length(roughPeakPos)
    %         plot(roughPeakPos(peak,1), roughPeakPos(peak,2), 'rx', 'MarkerSize', 15, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'none'); 
    %     end
    % end

    % Detect and mark fine peaks in the current chunk
    fineInRange     = (finePeakList(:,1) >= timeStart) & (finePeakList(:,1) <= timeEnd);
    if any(fineInRange)
        finePeakPos = finePeakList(fineInRange,:);
        for peak = 1:length(finePeakPos)
            plot(finePeakPos(peak,1)*pm.newTempRes, finePeakPos(peak,2), 'o', 'MarkerSize', 13, 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'none'); 
        end  
    end

    % Detect and mark breathStarts in the current chunk
    breathInRange = (breathStarts(1,:) >= timeStart) & (breathStarts(1,:) <= timeEnd);

    if any(breathInRange)
        breathPos = breathStarts(:,breathInRange);
        % Find the next peak (outside the window) to animate green area
        % between peaks.
        nextWindowPos = find(breathInRange,1,'last') + 1; 
        if nextWindowPos > size(breathStarts,2)
            return; % If we reach the end, return from function
        end
        breathPos(:,end+1) = breathStarts(:,nextWindowPos);

        for breath = 1:length(breathPos)-1
            % Plot the dot
            plot(breathPos(1,breath)*pm.newTempRes, breathPos(2,breath), 'o', 'MarkerSize', 13, 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'none');

            % Plot green area
            breathLength    = breathPos(1,breath+1) - breathPos(1,breath);
            breathStart     = breathPos(1,breath) + floor(breathLength*pm.breathRange(1));
            breathEnd       = breathPos(1,breath) + floor(breathLength*pm.breathRange(2));
            curHeight       = breathPos(2,breath);
            breathPoints    = breathEnd-breathStart+1;
            flatLine        = linspace(curHeight,curHeight,breathPoints);
            area((breathStart:breathEnd)*pm.newTempRes,flatLine,'FaceColor', 'green', 'EdgeColor', 'none', 'FaceAlpha', 0.2)
        end
    end

    ylim(ylimits);
    xlim([timeStart*pm.newTempRes,timeEnd*pm.newTempRes]);

    if step*stepSize < timepointsPerWindow*3 % Record specificed timerange into GIF
        GIFmaker(figHandle, gifFileName, 0.1);
    end
    hold off;
    pause(0.05);
end
end % End of function