%% CE653 - MONTE CARLO SIMULATION
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

%% 2. Choose the simulation size and reproducible random seed
% Keep the original 20 million draws for reproducible precision.
% Reduce totalSamples for a quick check; assess the reported confidence interval.
rng(653, 'twister');
totalSamples = 20000000;
batchSize = 100000;
assert(totalSamples >= 1 && totalSamples == floor(totalSamples) && batchSize >= 1);
failureCount = 0;
completedSamples = 0;
nonYieldCount = 0;
invalidCount = 0;
negativeLiveLoadCount = 0;
nonYieldFailureCount = 0;
meanMargin = 0;
sumSquaredDeviations = 0;
history = zeros(ceil(totalSamples / batchSize), 4);
runTimer = tic;
fprintf('\nMONTE CARLO: N = %d, batch size = %d\n', totalSamples, batchSize);
%% 3. Generate samples, evaluate g and count failures
for batchIndex = 1:ceil(totalSamples / batchSize)
    samplesInBatch = min(batchSize, totalSamples - completedSamples);
    samples = transform_uniform_samples(rand(samplesInBatch, numel(parameters.meanValues)), ...
                                        parameters);
    [margin, ~, ~, sampleDiagnostics] = evaluate_limit_state(samples, parameters);
    assert(all(isfinite(margin)), 'Nonfinite g encountered. Review the input model.');
    failureCount = failureCount + sum(margin <= 0);
    nonYieldCount = nonYieldCount + sum(sampleDiagnostics.nonYield);
    invalidCount = invalidCount + sum(sampleDiagnostics.invalid);
    negativeLiveLoadCount = negativeLiveLoadCount + sum(sampleDiagnostics.negativeQ);
    nonYieldFailureCount = nonYieldFailureCount + sum((margin <= 0) & sampleDiagnostics.nonYield);
    % Stable pooled sample moments, independent of memory batch size.
    % Pool batch moments without retaining all safety margins in memory.
    % sumSquaredDeviations is the accumulated sum about the running mean.
    batchMean = mean(margin);
    meanShift = batchMean - meanMargin;
    updatedSampleCount = completedSamples + samplesInBatch;
    sumSquaredDeviations = sumSquaredDeviations + sum((margin - batchMean).^2) + meanShift^2 * ...
        completedSamples * samplesInBatch / updatedSampleCount;
    meanMargin = meanMargin + meanShift * samplesInBatch / updatedSampleCount;
    completedSamples = updatedSampleCount;
    runningFailureProbability = failureCount / completedSamples;
    history(batchIndex, :) = [completedSamples, failureCount, runningFailureProbability, ...
                              reliability_index_from_probability(runningFailureProbability)];
    if mod(batchIndex, 20) == 0 || completedSamples == totalSamples
        fprintf('N=%9d failures=%5d Pf=%.6g beta=%.5f\n', completedSamples, failureCount, ...
                runningFailureProbability, ...
                reliability_index_from_probability(runningFailureProbability));
    end
end
%% 4. Estimate failure probability, beta and sampling uncertainty
failureProbability = failureCount / totalSamples;
beta = reliability_index_from_probability(failureProbability);
% Exact binomial Clopper-Pearson interval; only valid for independent MC.
[probabilityLowerBound, probabilityUpperBound] = binomial_interval(failureCount, totalSamples, ...
                                                                   0.05);
betaInterval = [reliability_index_from_probability(probabilityUpperBound), ...
                reliability_index_from_probability(probabilityLowerBound)];
standardError = sqrt(failureProbability * (1 - failureProbability) / totalSamples);
estimatorCOV = NaN;
if failureProbability > 0
    estimatorCOV = standardError / failureProbability;
end
%% 5. Store results and print the confidence intervals
results = struct('method', 'Crude Monte Carlo', 'parameters', parameters, 'N', totalSamples, ...
                 'failures', failureCount, 'Pf', failureProbability, 'beta', beta, 'Pf_CI95', ...
                 [probabilityLowerBound, probabilityUpperBound], ...
                 'beta_CI95', betaInterval, 'SE_Pf', standardError, 'COV_Pf_estimator', ...
                 estimatorCOV, ...
                 'sample_mean_g', meanMargin, 'sample_sd_g', sqrt(sumSquaredDeviations / ...
                                                                max(totalSamples - 1, 1)), ...
                 'history', history, 'non_yield_samples', nonYieldCount, 'invalid_samples', ...
                 invalidCount, ...
                 'failure_and_non_yield_samples', nonYieldFailureCount, 'elapsed_seconds', ...
                 toc(runTimer));
fprintf('\nMONTE CARLO RESULTS\n');
fprintf('Pf = %.8g; beta = %.8f\n', failureProbability, beta);
fprintf('Exact binomial 95%% CI(Pf): [%.8g, %.8g]\n', probabilityLowerBound, probabilityUpperBound);
fprintf('Transformed 95%% CI(beta): [%.6f, %.6f]\n', betaInterval);
fprintf('Estimated sampling COV(Pf) = %.2f%%\n', 100 * estimatorCOV);
fprintf('Sample mean(g)=%.6f; SD(g)=%.6f kN m\n', meanMargin, results.sample_sd_g);
fprintf('Failures outside steel-yield domain: %d\n', nonYieldFailureCount);
if failureCount == 0
    oneSidedUpperBound = 1 - 0.05^(1 / totalSamples);
    fprintf('ZERO failures: infinite point beta does NOT prove safety.\n');
    fprintf('One-sided 95%% upper Pf = %.8g; corresponding beta lower bound = %.5f\n', ...
            oneSidedUpperBound, reliability_index_from_probability(oneSidedUpperBound));
end
if failureCount < 100
    fprintf('Few failures: increase N for a more precise beta/Pf estimate (100 failures ~10%% COV).\n');
end
print_model_diagnostics(nonYieldCount, invalidCount, negativeLiveLoadCount, totalSamples);
fprintf('Elapsed: %.1f s. Output structure: results. No files written.\n', results.elapsed_seconds);

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

function samples = transform_uniform_samples(probabilities, parameters)
    % TRANSFORM_UNIFORM_SAMPLES  Inverse-CDF sampling from independent U(0,1).
    % First construct all normal variables; overwrite concrete and live load
    % with their lognormal and Gumbel inverse CDFs.

    % Protect log(0) and inverse-normal endpoints at machine precision only.
    probabilities = min(max(probabilities, realmin), 1 - eps);
    normalScores = -sqrt(2) * erfcinv(2 * probabilities);
    samples = bsxfun(@plus, parameters.meanValues, ...
                     bsxfun(@times, normalScores, parameters.standardDeviations));
    samples(:, 1) = exp(parameters.concreteLogMean ...
                        + parameters.concreteLogSD * normalScores(:, 1));
    samples(:, 10) = parameters.liveLoadLocation ...
        - parameters.liveLoadScale * log(-log(probabilities(:, 10)));
end

function beta = reliability_index_from_probability(failureProbability)
    % beta = -Phi^-1(Pf). erfcinv is built into MATLAB (no toolbox needed).
    % Pf=0 gives +Inf as a point estimate; it does not establish zero risk.
    beta = sqrt(2) * erfcinv(2 * failureProbability);
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

function [lowerBound, upperBound] = binomial_interval(failureCount, sampleCount, significance)
    % Exact two-sided Clopper-Pearson interval for independent Bernoulli draws.
    % The beta inverse CDF is available as base-MATLAB betaincinv.
    if failureCount == 0
        lowerBound = 0;
    else
        lowerBound = betaincinv(significance / 2, failureCount, sampleCount - failureCount + 1);
    end
    if failureCount == sampleCount
        upperBound = 1;
    else
        upperBound = betaincinv(1 - significance / 2, failureCount + 1, sampleCount - failureCount);
    end
end
