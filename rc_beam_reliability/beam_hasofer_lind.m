%% CE653 - HASOFER-LIND / FIRST-ORDER RELIABILITY METHOD (FORM)
% Design: 6 m span; 250 x 650 mm section; M25/Fe500; two 28 mm bars.
% Nominal line loads: G = 5 kN/m (excluding self-weight), Q = 30 kN/m.
% Run the COMPLETE script in MATLAB R2016b or later. No toolboxes required.
% Read the numbered sections first, then the local functions at the end.
% All basic inputs are independent. Both main bars share one diameter draw.
% The adopted yielded-steel capacity model and diagnostics are unchanged.
% Results are retained in the workspace structure named results.
% Edit define_beam_parameters() consistently in all four standalone files.

%% 1. Load the beam geometry and statistical input model
clear;
clc;
parameters = define_beam_parameters();
print_parameter_table(parameters);

%% 2. Set up the FORM problem in standard-normal coordinates
% Exact marginal transformation X_i = F_i^{-1}(Phi(U_i)).
% Solve min ||u|| subject to h(u)=g(T(u))=0 using damped HLRF.
% A merit-function line search prevents undamped HLRF oscillation.
% Multiple starts check for competing local design points (not a global proof).
marginScale = max(abs(evaluate_limit_state(parameters.meanValues, parameters)), 1);
scaledLimitState = ...
    @(normalCoordinates) evaluate_limit_state(transform_normal_coordinates(normalCoordinates, ...
                                                                           parameters), parameters) / marginScale;
numberOfVariables = numel(parameters.meanValues);
% Starts probe the origin, high live load, weak concrete, low resistance
% factor and weak steel. They check local convergence, not global optimality.
startingPoints = zeros(5, numberOfVariables);
startingPoints(2, 10) = 3;
startingPoints(3, 1) = -3;
startingPoints(4, 11) = -3;
startingPoints(5, 2) = -3;
options.maxIterations = 300;
options.tolH = 1e-8;
options.tolKKT = 1e-6;
options.fdStep = 1e-4;
solutions = cell(size(startingPoints, 1), 1);
distances = inf(size(startingPoints, 1), 1);
%% 3. Solve from each starting point and select the nearest converged point
fprintf('\nFORM MULTI-START CHECK\n');
for startIndex = 1:size(startingPoints, 1)
    solution = hlrf_solve(scaledLimitState, startingPoints(startIndex, :), options);
    solutions{startIndex} = solution;
    if solution.converged
        distances(startIndex) = norm(solution.u);
    end
    fprintf('Start %d: converged=%d, iterations=%d, ||u||=%.8f, |h|=%.3g\n', ...
            startIndex, solution.converged, solution.iterations, norm(solution.u), ...
            abs(scaledLimitState(solution.u)));
end
[designPointDistance, bestStartIndex] = min(distances);
if ~isfinite(designPointDistance)
    error('FORM did not converge. Inspect solutions and adjust the solver; no beta reported.');
end
solution = solutions{bestStartIndex};
designPointNormal = solution.u;
designPointPhysical = transform_normal_coordinates(designPointNormal, parameters);
[marginAtDesignPoint, resistanceAtDesignPoint, demandAtDesignPoint, designPointDiagnostics] = ...
    evaluate_limit_state(designPointPhysical, parameters);
%% 4. Check the design point and calculate beta, Pf and alpha
% The gradient points toward increasing g (the safer side).
% At a safe origin: u_star approximately equals -beta*alpha.
coarseGradient = calculate_central_gradient(scaledLimitState, designPointNormal, ...
                                            options.fdStep * ones(1, numberOfVariables));
designPointGradient = calculate_central_gradient(scaledLimitState, designPointNormal, 0.5 * ...
                                                 options.fdStep * ones(1, numberOfVariables));
relativeGradientDifference = norm(coarseGradient - designPointGradient) / ...
    max(norm(designPointGradient), eps);
assert(relativeGradientDifference < 1e-4, 'FORM gradient did not stabilize.');
alpha = designPointGradient / norm(designPointGradient); % points toward increasing safety
safeOriginSign = sign(scaledLimitState(zeros(1, numberOfVariables)));
if safeOriginSign == 0
    safeOriginSign = 1;
end
beta = safeOriginSign * designPointDistance;
failureProbability = 0.5 * erfc(beta / sqrt(2));
stationarityResidual = norm(designPointNormal - dot(designPointNormal, alpha) * alpha);
assert(abs(scaledLimitState(designPointNormal)) < options.tolH && stationarityResidual < 1e-5, ...
       'Final limit-state or stationarity residual is too large.');
%% 5. Store and print results
results = struct('method', 'Hasofer-Lind FORM', 'parameters', parameters, 'beta', beta, ...
                 'Pf_FORM', failureProbability, ...
                 'u_star', designPointNormal, 'x_star', designPointPhysical, 'g_star', ...
                 marginAtDesignPoint, 'R_star', resistanceAtDesignPoint, 'S_star', ...
                 demandAtDesignPoint, ...
                 'alpha_safe_direction', alpha, 'importance_alpha_squared', alpha.^2, ...
                 'stationarity_residual', stationarityResidual, 'gradient_check_relative', ...
                 relativeGradientDifference, ...
                 'history', solution.history, 'all_solutions', {solutions});
fprintf('\nHASOFER-LIND / FORM RESULTS\n');
fprintf('beta_FORM = %.8f; Pf_FORM = %.8g\n', beta, failureProbability);
fprintf('g(x*) = %.6g kN m; KKT residual = %.3g\n', marginAtDesignPoint, stationarityResidual);
fprintf('R(x*) = %.6f; S(x*) = %.6f kN m\n', resistanceAtDesignPoint, demandAtDesignPoint);
fprintf('\n%-11s %12s %12s %12s %12s\n', 'Variable', 'x*', 'u*', 'alpha', 'alpha^2');
for variableIndex = 1:numberOfVariables
    fprintf('%-11s %12.6g %12.6g %12.6g %12.6g\n', parameters.names{variableIndex}, ...
            designPointPhysical(variableIndex), designPointNormal(variableIndex), ...
            alpha(variableIndex), alpha(variableIndex)^2);
end
fprintf('\nalpha points toward safety; for this safe origin u* ~= -beta*alpha.\n');
fprintf('History columns: iteration, ||u||, h(u), HLRF residual.\n');
disp(solution.history);
print_model_diagnostics(sum(designPointDiagnostics.nonYield), ...
                        sum(designPointDiagnostics.invalid), sum(designPointDiagnostics.negativeQ), 1);
fprintf('FORM is a tangent-plane probability approximation, not exact integration.\n');
fprintf('The output structure is named results. No files are written.\n');

%% Local functions: input model, mechanics and numerical helpers

function parameters = define_beam_parameters()
    % DEFINE_BEAM_PARAMETERS  Agreed design and probability models.
    % Edit the nominal design here, then copy this function to all four scripts.
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
    nominalWidth = 250;
    nominalDepth = 650;
    nominalMainBarDiameter = 28;
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

function print_parameter_table(parameters)
    % Display arithmetic moments so the user can verify the actual run inputs.
    fprintf('\nIndependent basic variables (arithmetic mean and SD):\n');
    fprintf('%-11s %-13s %12s %12s %-8s\n', ...
            'Variable', 'Distribution', 'Mean', 'SD', 'Unit');
    for variableIndex = 1:numel(parameters.meanValues)
        fprintf('%-11s %-13s %12.6g %12.6g %-8s\n', ...
                parameters.names{variableIndex}, parameters.distributions{variableIndex}, ...
                parameters.meanValues(variableIndex), ...
                parameters.standardDeviations(variableIndex), parameters.units{variableIndex});
    end
    fprintf('Span = %.3f m; number of tension bars = %d.\n', ...
            parameters.span_m, parameters.numberOfBars);
    fprintf('Model: yielded-steel flexure with independent basic inputs.\n');
end

function print_model_diagnostics(nonYieldCount, invalidCount, negativeLiveLoadCount, sampleCount)
    % Report how often the adopted capacity model is used outside its domain.
    if invalidCount > 0
        warning('Beam:InvalidSamples', ...
                '%d / %d samples have invalid physical inputs. Samples retained.', ...
                invalidCount, sampleCount);
    end
    if negativeLiveLoadCount > 0
        warning('Beam:NegativeLiveLoad', ...
                '%d negative Gumbel live-load samples retained.', negativeLiveLoadCount);
    end
    if nonYieldCount > 0
        warning('Beam:YieldAssumption', ...
                ['%d / %d samples (%.5g%%) exceed the steel-yield strain limit. ', ...
                 'Results still use the adopted yielded-steel equation.'], ...
                nonYieldCount, sampleCount, 100 * nonYieldCount / sampleCount);
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
