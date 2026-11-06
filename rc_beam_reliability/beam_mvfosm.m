%% CE653 - MEAN-VALUE FIRST-ORDER SECOND-MOMENT (MVFOSM)
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

%% 2. Linearise the safety margin at the arithmetic mean
% mu_g ~= g(mu_X), Var(g) ~= grad_g * diag(sd.^2) * grad_g'.
% This is MVFOSM, NOT FORM, and does not use distribution shapes beyond moments.
limitStateFunction = @(x) evaluate_limit_state(x, parameters);
[marginAtMean, resistanceAtMean, demandAtMean, meanPointDiagnostics] = ...
    evaluate_limit_state(parameters.meanValues, parameters);
% Central differences perturb one physical input at a time.
% Halving the step checks numerical stability of the derivatives.
differenceSteps = 1e-4 * parameters.standardDeviations;
coarseGradient = calculate_central_gradient(limitStateFunction, parameters.meanValues, ...
                                            differenceSteps);
gradientAtMean = calculate_central_gradient(limitStateFunction, parameters.meanValues, ...
                                            differenceSteps / 2);
relativeGradientDifference = norm((coarseGradient - ...
                                   gradientAtMean) .* parameters.standardDeviations) / ...
    max(norm(gradientAtMean .* parameters.standardDeviations), eps);
assert(relativeGradientDifference < 1e-5, 'Finite-difference gradient did not stabilize.');
%% 3. Combine independent variance contributions and calculate beta
% For independent inputs, covariance cross-terms are zero.
varianceContributions = (gradientAtMean .* parameters.standardDeviations).^2;
marginVariance = sum(varianceContributions);
assert(marginVariance > 0, 'Zero first-order variance.');
marginSD = sqrt(marginVariance);
beta = marginAtMean / marginSD;
normalApproxFailureProbability = 0.5 * erfc(beta / sqrt(2));
varianceShares = varianceContributions / marginVariance;
%% 4. Store and print the approximation and sensitivity results
results = struct('method', 'MVFOSM', 'parameters', parameters, 'g_at_mean', marginAtMean, ...
                 'R_at_mean', resistanceAtMean, 'S_at_mean', demandAtMean, 'sigma_g_first_order', ...
                 marginSD, ...
                 'beta', beta, 'Pf_normal_approx', normalApproxFailureProbability, 'gradient_X', ...
                 gradientAtMean, ...
                 'variance_importance', varianceShares, 'gradient_check_relative', ...
                 relativeGradientDifference);
fprintf('\nMVFOSM RESULTS\n');
fprintf('R(mu) = %.6f; S(mu) = %.6f; g(mu) = %.6f kN m\n', resistanceAtMean, demandAtMean, ...
        marginAtMean);
fprintf('First-order SD(g) = %.6f kN m\n', marginSD);
fprintf('beta_MVFOSM = %.8f\n', beta);
fprintf('Phi(-beta_MVFOSM) = %.8g (normal-margin approximation only)\n', ...
        normalApproxFailureProbability);
fprintf('\n%-11s %15s %15s\n', 'Variable', 'dg/dX', 'Variance share');
for variableIndex = 1:numel(parameters.meanValues)
    fprintf('%-11s %15.6g %14.4f%%\n', parameters.names{variableIndex}, ...
            gradientAtMean(variableIndex), 100 * varianceShares(variableIndex));
end
print_model_diagnostics(sum(meanPointDiagnostics.nonYield), sum(meanPointDiagnostics.invalid), ...
                        sum(meanPointDiagnostics.negativeQ), 1);
fprintf('\nThe output structure is named results. No files are written.\n');

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
