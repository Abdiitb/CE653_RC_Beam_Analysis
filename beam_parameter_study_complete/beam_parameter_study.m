%% CE653 - PARAMETER STUDY OF FLEXURAL FORM RELIABILITY
% Vary b, overall D or bottom-bar diameter individually, holding the other
% nominal design parameters at the baseline. Both bottom bars share a draw.
% Probability-model rules remain unchanged; dimensional means/SDs and
% diameter SD (1% COV) update with their corresponding nominal dimensions.
% The original yielded-steel resistance equation is retained for every case.
% Two plots: separate panels in actual mm, and overlaid percentage changes.
% No toolboxes required. Run the COMPLETE script, MATLAB R2016b or later.
% The numeric figures supplied with this file were computed in Python;
% MATLAB will recompute using the same equations and a damped HLRF solver.

clear;
clc;

%% 1. Baseline and study ranges
baselineWidth = 250;
baselineDepth = 650;
baselineDiameter = 28;
targetBeta = 4.2;
widthValues = 200:25:450;
depthValues = 550:25:950;
diameterValues = 16:1:32; % Continuous sensitivity sweep, not a bar-size schedule.
percentageChanges = -20:5:20;
outputDirectory = fullfile(pwd, 'beam_parameter_study_outputs');
if ~exist(outputDirectory, 'dir')
    mkdir(outputDirectory);
end

%% 2. Analyse all three one-at-a-time sweeps
widthTable = run_sweep('Width', widthValues, baselineWidth, baselineDepth, baselineDiameter);
depthTable = run_sweep('Overall depth', depthValues, baselineWidth, baselineDepth, baselineDiameter);
diameterTable = run_sweep('Bar diameter', diameterValues, baselineWidth, baselineDepth, baselineDiameter);
resultsTable = [widthTable; depthTable; diameterTable];
writetable(resultsTable, fullfile(outputDirectory, 'beam_parameter_study_data.csv'));

%% 3. Plot separate panels using actual nominal dimensions (mm)
figurePanels = figure('Color', 'w', 'Position', [100, 100, 1400, 550]);
colors = [0.145, 0.388, 0.651; 0.086, 0.510, 0.451; 0.651, 0.290, 0.467];
tables = {widthTable, depthTable, diameterTable};
labels = {'Width b (mm)', 'Overall depth D (mm)', 'Bottom-bar diameter (mm)'};
fixedInputs = {'D=650 mm; two 28 mm bars', 'b=250 mm; two 28 mm bars', ...
               'b=250 mm; D=650 mm; two bars'};
baselineParameters = [baselineWidth, baselineDepth, baselineDiameter];
baselineAnalysis = analyse_form(baselineWidth, baselineDepth, baselineDiameter);
for parameterIndex = 1:3
    subplot(1, 3, parameterIndex);
    currentTable = tables{parameterIndex};
    plot(currentTable.Value_mm, currentTable.Beta_FORM, '-o', ...
         'Color', colors(parameterIndex, :), 'LineWidth', 1.8);
    hold on;
    plot([min(currentTable.Value_mm), max(currentTable.Value_mm)], ...
         [targetBeta, targetBeta], '--', 'Color', [0.75, 0.27, 0.2]);
    plot(baselineParameters(parameterIndex), baselineAnalysis.beta, ...
         'kp', 'MarkerFaceColor', 'k', 'MarkerSize', 11);
    notAcceptable = ~currentTable.DuctilityPass | ~currentTable.CoverPass;
    plot(currentTable.Value_mm(notAcceptable), currentTable.Beta_FORM(notAcceptable), ...
         'kx', 'MarkerSize', 9, 'LineWidth', 1.5);
    xlabel(labels{parameterIndex});
    ylabel('FORM reliability index, beta');
    title(fixedInputs{parameterIndex});
    grid on;
    box off;
end
print(figurePanels, fullfile(outputDirectory, 'beam_beta_comparison.png'), '-dpng', '-r300');
set(figurePanels, 'PaperPositionMode', 'auto');
print(figurePanels, fullfile(outputDirectory, 'beam_beta_comparison.pdf'), '-dpdf', '-bestfit');

%% 4. Recompute at equal percentage changes for a fair overlaid comparison
% Different dimension magnitudes cannot be compared fairly on one raw-mm
% axis. At +10%, for example: b=275 OR D=715 OR diameter=30.8 mm.
widthOverlay = run_sweep('Width', baselineWidth * (1 + percentageChanges / 100), ...
                         baselineWidth, baselineDepth, baselineDiameter);
depthOverlay = run_sweep('Overall depth', baselineDepth * (1 + percentageChanges / 100), ...
                         baselineWidth, baselineDepth, baselineDiameter);
diameterOverlay = run_sweep('Bar diameter', baselineDiameter * (1 + percentageChanges / 100), ...
                            baselineWidth, baselineDepth, baselineDiameter);
overlayTables = {widthOverlay, depthOverlay, diameterOverlay};
overlayData = [widthOverlay; depthOverlay; diameterOverlay];
overlayData.PercentageChange = repmat(percentageChanges(:), 3, 1);
writetable(overlayData, fullfile(outputDirectory, 'beam_parameter_overlay_data.csv'));

figureOverlay = figure('Color', 'w', 'Position', [200, 100, 850, 600]);
hold on;
% Use a cell array of handles for compatibility with MATLAB graphics objects.
curveHandles = cell(1, 3);
for parameterIndex = 1:3
    currentTable = overlayTables{parameterIndex};
    curveHandles{parameterIndex} = plot(percentageChanges, currentTable.Beta_FORM, '-o', ...
                                        'Color', colors(parameterIndex, :), 'LineWidth', 1.8);
    notAcceptable = ~currentTable.DuctilityPass | ~currentTable.CoverPass;
    plot(percentageChanges(notAcceptable), currentTable.Beta_FORM(notAcceptable), ...
         'kx', 'MarkerSize', 10, 'LineWidth', 1.5);
end
targetHandle = plot([min(percentageChanges), max(percentageChanges)], ...
                    [targetBeta, targetBeta], '--', 'Color', [0.75, 0.27, 0.2], 'LineWidth', 1.4);
plot(0, baselineAnalysis.beta, 'kp', 'MarkerFaceColor', 'k', 'MarkerSize', 12);
legend([curveHandles{1}, curveHandles{2}, curveHandles{3}, targetHandle], ...
       {'Width (baseline 250 mm)', 'Overall depth (baseline 650 mm)', ...
        'Bar diameter (baseline 28 mm)', 'Target beta = 4.2'}, 'Location', 'northwest');
xlabel('Change from corresponding nominal baseline (%)');
ylabel('FORM reliability index, beta');
title({'One-at-a-time parameter comparison', ...
       'Crosses: nominal ductility or fixed-cover constraint fails'});
grid on;
box off;
print(figureOverlay, fullfile(outputDirectory, 'beam_beta_overlay.png'), '-dpng', '-r300');
set(figureOverlay, 'PaperPositionMode', 'auto');
print(figureOverlay, fullfile(outputDirectory, 'beam_beta_overlay.pdf'), '-dpdf', '-bestfit');

%% 5. Keep data and settings in the workspace
results.baselineBeta = baselineAnalysis.beta;
results.parameterSweeps = resultsTable;
results.normalizedComparison = overlayData;
results.targetBeta = targetBeta;
results.baselineGeometry = baselineParameters;
fprintf('\nBaseline beta = %.6f. Figures and CSVs saved in:\n%s\n', ...
        baselineAnalysis.beta, outputDirectory);
fprintf('These are model-based sensitivity curves, not a schedule of approved designs.\n');
fprintf('A negative signed beta indicates failure at the standard-normal origin.\n');
fprintf('Only the selected nominal parameter changes in each curve.\n');

%% Local functions

function resultTable = run_sweep(parameterName, values, baseWidth, baseDepth, baseDiameter)
    pointCount = numel(values);
    betaValues = zeros(pointCount, 1);
    probabilityValues = zeros(pointCount, 1);
    geometry = zeros(pointCount, 3);
    ductilityPass = false(pointCount, 1);
    coverPass = false(pointCount, 1);
    designPointYieldPass = false(pointCount, 1);
    for pointIndex = 1:pointCount
        width = baseWidth;
        depth = baseDepth;
        diameter = baseDiameter;
        switch parameterName
            case 'Width'
                width = values(pointIndex);
            case 'Overall depth'
                depth = values(pointIndex);
            case 'Bar diameter'
                diameter = values(pointIndex);
        end
        analysis = analyse_form(width, depth, diameter);
        betaValues(pointIndex) = analysis.beta;
        probabilityValues(pointIndex) = analysis.Pf;
        geometry(pointIndex, :) = [width, depth, diameter];
        area = 2 * pi * diameter^2 / 4;
        effectiveDepth = depth - 30 - 8 - diameter / 2;
        neutralAxis = 0.87 * 500 * area / (0.36 * 25 * width);
        ductilityPass(pointIndex) = neutralAxis / effectiveDepth <= 0.46;
        coverPass(pointIndex) = diameter <= 30;
        designPointYieldPass(pointIndex) = ~analysis.diagnostics.nonYield;
        fprintf('%s = %.3f mm: beta = %.6f\n', parameterName, values(pointIndex), analysis.beta);
    end
    resultTable = table(repmat({parameterName}, pointCount, 1), values(:), ...
                        geometry(:, 1), geometry(:, 2), geometry(:, 3), betaValues, probabilityValues, ...
                        ductilityPass, coverPass, designPointYieldPass, ...
                        'VariableNames', {'VariedParameter', 'Value_mm', 'Width_mm', 'Depth_mm', ...
                          'Diameter_mm', 'Beta_FORM', 'Pf_FORM', 'DuctilityPass', 'CoverPass', 'DesignPointYieldPass'});
end

function margin = margin_at_normal_point(normalPoint, nominalWidth, nominalDepth, nominalMainBarDiameter)
    parameters = define_beam_parameters(nominalWidth, nominalDepth, nominalMainBarDiameter);
    physicalPoint = transform_normal_coordinates(normalPoint, parameters);
    margin = evaluate_limit_state(physicalPoint, parameters);
end

function analysis = analyse_form(nominalWidth, nominalDepth, nominalMainBarDiameter)
    % Find minimum ||u|| subject to g(T(u))=0, independently of target design.
    parameters = define_beam_parameters(nominalWidth, nominalDepth, nominalMainBarDiameter);
    numberOfVariables = numel(parameters.meanValues);
    marginScale = max(abs(evaluate_limit_state(parameters.meanValues, parameters)), 1);
    scaledFunction = @(u) margin_at_normal_point(u, nominalWidth, nominalDepth, nominalMainBarDiameter) / marginScale;
    startingPoints = zeros(5, numberOfVariables);
    startingPoints(2, 10) = 3;
    startingPoints(3, 1) = -3;
    startingPoints(4, 11) = -3;
    startingPoints(5, 2) = -3;
    options.maxIterations = 300;
    options.tolH = 1e-8;
    options.tolKKT = 1e-6;
    options.fdStep = 1e-4;
    bestDistance = Inf;
    bestSolution = [];
    for startIndex = 1:size(startingPoints, 1)
        solution = hlrf_solve(scaledFunction, startingPoints(startIndex, :), options);
        if solution.converged && norm(solution.u) < bestDistance
            bestDistance = norm(solution.u);
            bestSolution = solution;
        end
    end
    assert(~isempty(bestSolution), 'Independent FORM analysis failed to converge.');
    normalPoint = bestSolution.u;
    physicalPoint = transform_normal_coordinates(normalPoint, parameters);
    [margin, resistance, demand, diagnostics] = evaluate_limit_state(physicalPoint, parameters);
    gradient = calculate_central_gradient(scaledFunction, normalPoint, ...
                                          options.fdStep * ones(1, numberOfVariables));
    failureDirection = -gradient / norm(gradient);
    stationarityResidual = norm(normalPoint - dot(normalPoint, failureDirection) * failureDirection);
    assert(abs(margin) < 1e-5 && stationarityResidual < 1e-5, ...
           'FORM point failed the boundary or stationarity check.');
    originSign = sign(scaledFunction(zeros(1, numberOfVariables)));
    analysis.beta = originSign * bestDistance;
    analysis.Pf = 0.5 * erfc(analysis.beta / sqrt(2));
    analysis.nominalDepth = nominalDepth;
    analysis.uStar = normalPoint;
    analysis.xStar = physicalPoint;
    analysis.margin = margin;
    analysis.resistance = resistance;
    analysis.demand = demand;
    analysis.failureDirection = failureDirection;
    analysis.stationarityResidual = stationarityResidual;
    analysis.diagnostics = diagnostics;
    analysis.history = bestSolution.history;
end

function parameters = define_beam_parameters(nominalWidth, nominalDepth, nominalMainBarDiameter)
    % DEFINE_BEAM_PARAMETERS  Agreed design and probability models.
    % All statistics are retained; only nominalDepth is varied by the caller.
    % meanValues and standardDeviations are ARITHMETIC moments, including for
    % lognormal concrete strength and Gumbel live load.
    %
    % Column order in every sample matrix:
    % 1 concrete strength; 2 yield strength; 3 width; 4 overall depth; 5 cover;
    % 6 main-bar diameter; 7 stirrup diameter; 8 RC unit weight;
    % 9 superimposed dead load; 10 live load; 11 resistance model factor.

    % Nominal design quantities: lengths in mm unless explicitly stated.
    parameters.span_m = 6;
    parameters.numberOfBars = 2;
    % nominalDepth is supplied by the design iteration, in mm.
    nominalStirrupDiameter = 8;
    nominalDeadLoad = 5;              % kN/m, EXCLUDING beam self-weight
    nominalLiveLoad = 30;             % kN/m
    characteristicConcreteStrength = 25; % MPa (M25)

    % These constants are used only to check the steel-yield assumption.
    parameters.steelModulus_MPa = 200000;
    parameters.ultimateConcreteStrain = 0.0035;

    % Concrete: lognormal, 15% COV, with 25 MPa at the fifth percentile.
    % Adopted transfer assumption; not a fitted M25 test dataset.
    concreteCOV = 0.15;
    parameters.concreteLogSD = sqrt(log(1 + concreteCOV^2));
    parameters.concreteLogMean = log(characteristicConcreteStrength) ...
        + 1.6448536269514722 * parameters.concreteLogSD;
    meanConcreteStrength = exp(parameters.concreteLogMean ...
                               + 0.5 * parameters.concreteLogSD^2);
    sdConcreteStrength = concreteCOV * meanConcreteStrength;

    % JCSS external-dimension model: mean = 1.003*nominal;
    % SD = 4 + 0.006*nominal, with dimensions in mm.
    meanWidth = 1.003 * nominalWidth;
    sdWidth = 4 + 0.006 * nominalWidth;
    meanDepth = 1.003 * nominalDepth;
    sdDepth = 4 + 0.006 * nominalDepth;

    % Normal diameter model with 1% COV. Both main bars share ONE diameter
    % draw, preserving the original common-diameter assumption.
    sdMainBarDiameter = 0.01 * nominalMainBarDiameter;
    sdStirrupDiameter = 0.01 * nominalStirrupDiameter;

    % Correia et al. (2022), Table 2: transferred statistical load model.
    meanDeadLoad = 1.05 * nominalDeadLoad;
    sdDeadLoad = 0.10 * meanDeadLoad;
    liveLoadCOV = 0.25;
    meanLiveLoad = nominalLiveLoad / (1 + 0.35 * liveLoadCOV);
    sdLiveLoad = liveLoadCOV * meanLiveLoad;

    parameters.names = {'C', 'F', 'b', 'D', 'c', 'phi_l', 'phi_st', ...
                        'gamma_RC', 'G', 'Q', 'theta_R'};
    parameters.units = {'MPa', 'MPa', 'mm', 'mm', 'mm', 'mm', 'mm', ...
                        'kN/m^3', 'kN/m', 'kN/m', '-'};

    % Assemble the 11 basic variables in the column order listed above.
    parameters.meanValues = [meanConcreteStrength, 560, meanWidth, ...
                             meanDepth, 30, nominalMainBarDiameter, nominalStirrupDiameter, ...
                             25, meanDeadLoad, meanLiveLoad, 1.01];
    parameters.standardDeviations = [sdConcreteStrength, 30, sdWidth, ...
                                     sdDepth, 3, sdMainBarDiameter, sdStirrupDiameter, ...
                                     1.25, sdDeadLoad, sdLiveLoad, 0.0465];
    parameters.distributions = {'Lognormal', 'Normal', 'Normal', 'Normal', ...
                                'Normal', 'Normal', 'Normal', 'Normal', 'Normal', 'Gumbel-max', ...
                                'Normal'};

    % Gumbel parameters are NOT the arithmetic mean and SD.
    % scale = SD*sqrt(6)/pi; location = mean - EulerConstant*scale.
    parameters.liveLoadScale = sdLiveLoad * sqrt(6) / pi;
    parameters.liveLoadLocation = meanLiveLoad ...
        - 0.5772156649015329 * parameters.liveLoadScale;
end

function [margin, resistance, demand, diagnostics] = evaluate_limit_state(samples, parameters)
    % EVALUATE_LIMIT_STATE  Flexural safety margin g = R - S, in kN*m.
    % samples: one row per realisation, 11 columns in the documented order.
    % A negative or zero margin counts as failure. No partial safety factors
    % (1.5 or 0.87) are included in this probabilistic model.

    concreteStrength = samples(:, 1);       % MPa
    steelYieldStrength = samples(:, 2);     % MPa
    width = samples(:, 3);                  % mm
    overallDepth = samples(:, 4);           % mm
    clearCover = samples(:, 5);             % mm, to outside of stirrup
    mainBarDiameter = samples(:, 6);        % mm
    stirrupDiameter = samples(:, 7);        % mm
    rcUnitWeight = samples(:, 8);           % kN/m^3, includes steel allowance
    deadLoad = samples(:, 9);               % kN/m, excludes beam self-weight
    liveLoad = samples(:, 10);              % kN/m
    resistanceFactor = samples(:, 11);      % dimensionless

    % Derived geometry: use the SAME sampled dimensions in every equation.
    steelArea = parameters.numberOfBars * pi * mainBarDiameter.^2 / 4;
    effectiveDepth = overallDepth - clearCover ...
        - stirrupDiameter - mainBarDiameter / 2;
    selfWeight = rcUnitWeight .* width .* overallDepth * 1e-6;

    % Adopted yielded-steel resistance equation. The 0.77 coefficient is
    % retained exactly from the original model. MPa*mm^3 = N*mm;
    % multiplying by 1e-6 converts N*mm to kN*m.
    steelTensionForce = steelYieldStrength .* steelArea;
    leverArm = effectiveDepth ...
        - 0.77 * steelTensionForce ./ (width .* concreteStrength);
    resistance = resistanceFactor .* steelTensionForce .* leverArm * 1e-6;

    % Simply supported beam under a uniformly distributed total load.
    totalLoad = deadLoad + liveLoad + selfWeight;
    demand = totalLoad * parameters.span_m^2 / 8;
    margin = resistance - demand;

    if nargout > 3
        % Diagnostics flag model-domain issues; they do NOT discard samples
        % or switch to another capacity equation.
        diagnostics.invalid = any(~isfinite(samples), 2) ...
            | any(samples(:, [1:9, 11]) <= 0, 2) | effectiveDepth <= 0;
        diagnostics.negativeQ = liveLoad < 0;
        neutralAxisDepth = steelTensionForce ...
            ./ (0.546 * concreteStrength .* width);
        diagnostics.yieldLimit = parameters.ultimateConcreteStrain ...
            ./ (parameters.ultimateConcreteStrain ...
                + steelYieldStrength / parameters.steelModulus_MPa);
        diagnostics.xuOverD = neutralAxisDepth ./ effectiveDepth;
        diagnostics.nonYield = diagnostics.xuOverD > diagnostics.yieldLimit;
        diagnostics.As = steelArea;
        diagnostics.d = effectiveDepth;
        diagnostics.wSelf = selfWeight;
    end
end

function samples = transform_normal_coordinates(normalCoordinates, parameters)
    % TRANSFORM_NORMAL_COORDINATES  FORM transformation X_i = F_i^-1(Phi(U_i)).
    % The origin U=0 maps to each variable's median, not necessarily its mean.

    samples = bsxfun(@plus, parameters.meanValues, ...
                     bsxfun(@times, normalCoordinates, parameters.standardDeviations));
    samples(:, 1) = exp(parameters.concreteLogMean ...
                        + parameters.concreteLogSD * normalCoordinates(:, 1));
    liveLoadProbability = 0.5 * erfc(-normalCoordinates(:, 10) / sqrt(2));
    liveLoadProbability = min(max(liveLoadProbability, realmin), 1 - eps);
    samples(:, 10) = parameters.liveLoadLocation ...
        - parameters.liveLoadScale * log(-log(liveLoadProbability));
end

function gradient = calculate_central_gradient(functionHandle, point, stepSizes)
    % CALCULATE_CENTRAL_GRADIENT  One numerical partial derivative per variable.
    % dg/dx_i ~= [g(x+h_i e_i) - g(x-h_i e_i)] / (2*h_i).

    gradient = zeros(size(point));
    for variableIndex = 1:numel(point)
        pointPlus = point;
        pointMinus = point;
        pointPlus(variableIndex) = point(variableIndex) + stepSizes(variableIndex);
        pointMinus(variableIndex) = point(variableIndex) - stepSizes(variableIndex);
        gradient(variableIndex) = (functionHandle(pointPlus) ...
                                   - functionHandle(pointMinus)) / (2 * stepSizes(variableIndex));
    end
end

function solution = hlrf_solve(scaledLimitState, currentPoint, options)
    % HLRF_SOLVE  Damped Hasofer-Lind-Rackwitz-Fiessler iteration.
    % At each point u, linearise h(v) ~= h(u) + grad(h)*(v-u).
    % The nearest point on that plane is
    % v = [(grad(h)*u - h(u))/||grad(h)||^2] * grad(h).
    % Move toward v, reducing the step until a merit function improves.
    % Merit = 0.5*||u||^2 + penalty*abs(h(u)); it balances distance and
    % distance from the failure surface. Both convergence checks must pass.
    % history columns: iteration, ||u||, h(u), ||v-u||.
    penaltyWeight = 1;
    iterationHistory = zeros(options.maxIterations + 1, 4);
    converged = false;
    reason = 'maximum iterations';
    for iteration = 0:options.maxIterations
        currentMargin = scaledLimitState(currentPoint);
        gradient = calculate_central_gradient(scaledLimitState, currentPoint, options.fdStep * ...
                                              ones(size(currentPoint)));
        squaredGradientNorm = dot(gradient, gradient);
        if ~isfinite(currentMargin) || ~isfinite(squaredGradientNorm) || squaredGradientNorm <= eps^2
            reason = 'nonfinite function or zero gradient';
            break
        end

        % Project the origin onto the current tangent plane.
        linearisedDesignPoint = ((dot(gradient, currentPoint) - currentMargin) / ...
                                 squaredGradientNorm) * gradient;
        searchDirection = linearisedDesignPoint - currentPoint;
        iterationHistory(iteration + 1, :) = [iteration, norm(currentPoint), currentMargin, ...
                                              norm(searchDirection)];
        if abs(currentMargin) < options.tolH && norm(searchDirection) < options.tolKKT
            converged = true;
            reason = 'converged';
            break
        end
        if iteration == options.maxIterations
            break
        end

        % Choose a penalty large enough to favour approaching h=0.
        penaltyWeight = max(penaltyWeight, 2 * norm(currentPoint) / sqrt(squaredGradientNorm) + 1);
        meritSlope = dot(currentPoint, searchDirection) - penaltyWeight * abs(currentMargin);
        if meritSlope >= 0 && abs(currentMargin) > 1e-14
            penaltyWeight = max(penaltyWeight, (dot(currentPoint, searchDirection) + ...
                                                norm(searchDirection)^2) / abs(currentMargin) + 1);
            meritSlope = dot(currentPoint, searchDirection) - penaltyWeight * abs(currentMargin);
        end
        currentMerit = 0.5 * dot(currentPoint, currentPoint) + penaltyWeight * abs(currentMargin);
        accepted = false;
        stepLength = 1;

        % Armijo backtracking: halve the proposed step until it improves merit.
        for backtrackIndex = 1:40
            trialPoint = currentPoint + stepLength * searchDirection;
            trialMargin = scaledLimitState(trialPoint);
            trialMerit = 0.5 * dot(trialPoint, trialPoint) + penaltyWeight * abs(trialMargin);
            if isfinite(trialMerit) && trialMerit <= currentMerit + 1e-4 * stepLength * meritSlope
                accepted = true;
                break
            end
            stepLength = stepLength / 2;
        end
        if ~accepted
            reason = 'line search failed';
            break
        end
        currentPoint = trialPoint;
    end
    solution = struct('u', currentPoint, 'converged', converged, 'iterations', iteration, ...
                      'reason', reason, ...
                      'history', iterationHistory(1:iteration + 1, :));
end
