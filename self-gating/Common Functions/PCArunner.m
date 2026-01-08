function [allCoeff, allScore, allExplained] = PCArunner(inputMidlines, pm, visualSwitch)
% Takes in an inputMatrix (e.g. midlines from cineReader.m [4x96x1152] [pixels, cardiac phase,
% coils]) and calculates PCA stuff.


disp('Section 2: PCA');
tic;

% PCA must be run on real and imag separately, so we split them
realMatrix  = real(inputMidlines);
imagMatrix  = imag(inputMidlines);

rawMidlines = []; % The vector to run through PCA
for coil = 1:pm.coilNum
    realMidlines        = squeeze(realMatrix(coil,:,:));
    imagMidlines        = squeeze(imagMatrix(coil,:,:));

    rawMidlines         = [rawMidlines, realMidlines', imagMidlines'];
end

standardizedMidlines    = zscore(rawMidlines,0,2);

% Run PCA
[allCoeff, allScore, allExplained] = pca(standardizedMidlines);
% 'coeff' contains the principal component coefficients (loadings)
% 'score' contains the principal component scores
% 'explained' contains the percentage of total variance explained by each PC

% Make variable first index and time second
allCoeff        = allCoeff';
allScore        = allScore';
allExplained    = allExplained';

% Analyze the first principal component
allScore1 = allScore(1, :); % The first PC score
allScore2 = allScore(2, :);
allScore3 = allScore(3, :);

if visualSwitch == true
    figNumStart = 6;

    % Visualize the explained variance
    figure(figNumStart);
    pareto(allExplained);
    xlabel('Principal Component');
    ylabel('Variance Explained (%)');
    title('All lines explained');
    
    % Visualize the first principal component
    allAxis = (1:size(inputMidlines,3));
    
    % Figures
    if length(allScore1) >= 250
        currentRange = 1:250;
    else
        currentRange = 1:length(allScore1);
    end
    figure(figNumStart+1);
    plot(allAxis(currentRange),allScore1(currentRange));
    hold on;
    xlabel('Time Instance');
    ylabel('Principal Component 1');
    title('All lines: First Principal Component Reflecting Cardiac Phase');
    legend('Location','best');
    hold off;
    
    if length(allScore1) >= 150
        currentRange = 1:150;
    else
        currentRange = 1:length(allScore1);
    end
    figure(figNumStart+2);
    plot(allAxis(currentRange),allScore2(currentRange));
    hold on;
    xlabel('Time Instance');
    ylabel('Principal Component 2');
    title('All lines: Second Principal Component Reflecting Something Else');
    legend('Location','best');
    hold off;
    
    figure(figNumStart+3);
    plot(allAxis(currentRange),allScore3(currentRange));
    hold on;
    xlabel('Time Instance');
    ylabel('Principal Component 3');
    title('All lines: Third Principal Component Reflecting Something Else');
    legend('Location','best');
    hold off;
    
    chosenCoeff = 1; % Based on the previous figures, select the component representing cardiac motion
    allOnAll    = standardizedMidlines * allCoeff(chosenCoeff,:)';
    
    figure(figNumStart+4);
    plot(allAxis,allOnAll);
    hold on;
    xlabel('Timepoint');
    ylabel('Magnitude');
    title([num2str(chosenCoeff) ' PC - Full data']);
    legend('Location','best');
    hold off;
end

disp(['Section 2: PCArunner - Finished in ' num2str(toc) ' seconds.' newline]);
end % of main function

function fixedList = signFlipper(pcList)
    % Flips the sign of PCs so that they point the same way in graphics.
    % Sign of PCs appear to be psuedo-random in any case.
    fixedList = pcList;
    for i = 1:length(pcList)
        tempPc = pcList{i};
        if tempPc(round(length(tempPc)*0.5)) < 0 % [FIKS] Manual choice of where to evaluate flipping-need
            tempPc = -tempPc;
            fixedList{i} = tempPc;
        end
    end
end % End of signFlipper