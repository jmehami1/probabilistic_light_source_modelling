function whiteCorr = WhiteReflectanceCorrection(S_estimated, S_whiteHypCube, figDir, methodName)

S_ratio = S_whiteHypCube./S_estimated;
S_ratio(isnan(S_ratio)) = 0;
S_ratio(isinf(S_ratio)) = 0;

% Reflectance Correction factor should be calculated per pixel (averaged over lines and bands)
whiteCorr = squeeze(median(S_ratio, [2,3]));

%fit polynomial curve to alpha correction to get correct reflectance values
[xData, yData] = prepareCurveData( [], whiteCorr );

% Set up fittype and options.
ft = fittype( 'poly8' );
excludedPoints = excludedata( xData, yData, 'Indices', find(yData <= 0) );
opts = fitoptions( 'Method', 'LinearLeastSquares' );
opts.Normalize = 'on';
opts.Robust = 'Bisquare';
opts.Exclude = excludedPoints;

% Fit model to data.
whiteCorrFit = fit( xData, yData, ft, opts );
whiteCorr = whiteCorrFit(1:length(whiteCorr));

figure();
plot(whiteCorrFit, xData, yData); hold on;
grid on;
ylabel('Reflectance Scale Correction Factor')
xlabel('Line-scan Pixels');
title('Reflectance scale correction factor across line-scan pixels');
savefig([figDir, filesep, 'reflectance_correction_', methodName, '.fig']);


end

