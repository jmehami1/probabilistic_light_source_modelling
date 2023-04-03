% function [mu, varMu, hypOpt] = LightSrcOptmGP(meanType, trainingData, testingX, sigmaNoise, downSamp, iter, verboseFlag, varargin)
% function [mu, varMu, hypOpt, trainingDownXY] = LightSrcOptmGP(meanType, trainingXY, testingX, sigmaNoise, varargin)


function [testingY, var_testingY, hypOptStruct] = LightSrcOptmGP(testingX, varargin)
% Builds the intensity field of a non-isotropic disk light source model
% given data using Gaussian Processes with different mean-function options.
%
% INPUTS:
%       testingX - Testing data array organised as columns [radius, theta]
%
%       ****Arugments for training and querying****
%       trainingXY - training data array organised as columns [radius, theta, radiant intensity magnitude]
%       meanType - Type of GP mean function
%                       0 ~ zero mean
%                       1 ~ constant mean
%                       2 ~ light source mean function (DEFAULT)
%                       3 ~ light source mean function with exponetial hyperparameters
%       sigmaNoise - STD of Gaussian noise in measurements
%       downSamp - Down-sampling rate applied to training data
%       iter - maximum number of iterations during optimisation
%       verboseFlag - verbose flag to show output to terminal (true/false)
%
%       ****Arguments for querying only****
%       hypOptStruct - struct containing optimised hyperparameter struct from prior training of GP
%           model
% 
% OUTPUTS:
%       testingY - Queried radiant intensity mean of testing data from optimised
%           GP model
%       var_testingY - Queried radiant intensity variance of testing data from
%           optimised GP model
%       hypOptStruct - struct of parameters that store optimised
%           hyperparameter struct and other GP information
% 
% Author: Jasprabhjit Mehami, 13446277

%only querying optimised model
if nargin == 2
    hypOptStruct = varargin{1};
    
    hypOpt = hypOptStruct.hypOpt;
    meanfunc = hypOptStruct.meanfunc;
    covfunc = hypOptStruct.covfunc;
    likfunc = hypOptStruct.likfunc;
    
    trainingDownXY = hypOptStruct.trainingDownXY;
    xGP = trainingDownXY(:,1:2);
    yGP = trainingDownXY(:,3);
    
    %training of GP and querying optimised model
elseif nargin == 7
    trainingXY = varargin{1};
    meanType = varargin{2};
    sigmaNoise = varargin{3};
    downSamp = varargin{4};
    iter = varargin{5};
    verboseFlag = varargin{6};
    
    covfunc = @covSEiso; %covariance function
    likfunc = @likGauss; %gaussian likelihood
    
    switch (meanType)
        %zero-mean
        case 0
            meanfunc = [];
            meanHyper = [];
            %constant mean
        case 1
            meanfunc = @meanConst;
            meanHyper = 1;
            %light source mean function
        case 2
            meanfunc = @meanLightSrc;
            meanHyper = [1,1,1];
        case 3
            meanfunc = @meanLightSrcExp;
            meanHyper = [1,1,1];
        otherwise
            meanfunc = @meanLightSrc;
            meanHyper = [1,1,1];
    end
    
    %initial hyperparameter struct
    hyp = struct('mean', meanHyper, 'cov', [0,0], 'lik', log(sigmaNoise));
    
    %downsample data
    trainingDownXY = downsample(trainingXY, downSamp);
    xGP = trainingDownXY(:,1:2);
    yGP = trainingDownXY(:,3);
    
    hypOpt = minimize(hyp, @gp, -iter, verboseFlag, @infGaussLik, meanfunc, covfunc, likfunc, xGP, yGP);
    
    %save all parameters into struct
    hypOptStruct.hypOpt = hypOpt;
    hypOptStruct.meanfunc = meanfunc;
    hypOptStruct.covfunc = covfunc;
    hypOptStruct.likfunc = likfunc;
    hypOptStruct.trainingDownXY = trainingDownXY;
else
    error('There should be either 2 or 7 input arguments for querying only or training and querying respectively.');
end

if nargout > 1
    %use the optimised hyperparameter struct to get the mean and variance of
    %given testing points
    [testingY, var_testingY] = gp(hypOpt, @infGaussLik, meanfunc, covfunc, likfunc, xGP, yGP, testingX);

else
    testingY = gp(hypOpt, @infGaussLik, meanfunc, covfunc, likfunc, xGP, yGP, testingX);
end

%no imaginery elements
if ~isreal(testingY)
    testingY = real(testingY);
end

end

%
% %only querying optimised model
% if nargin == 5
%     xGP = trainingXY(:,1:2);
%     yGP = trainingXY(:,3);
%
%     hypOpt = varargin{1};
%     %training of GP and querying optimised model
% elseif nargin == 7
%     %hyperparameter struct
%     hyp = struct('mean', meanHyper, 'cov', [0,0], 'lik', log(sigmaNoise));
%
%     downSamp = varargin{1};
%     iter = varargin{2};
%     verboseFlag = varargin{3};
%
%     %downsample data
%     trainingDownXY = downsample(trainingXY, downSamp);
%     xGP = trainingDownXY(:,1:2);
%     yGP = trainingDownXY(:,3);
%
%     hypOpt = minimize(hyp, @gp, -iter, verboseFlag, @infGaussLik, meanfunc, covfunc, likfunc, xGP, yGP);
% else
%     error('There should be either 5 or 7 input arguments');
% end
