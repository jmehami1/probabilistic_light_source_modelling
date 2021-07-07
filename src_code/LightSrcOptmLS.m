function [optPhi,resNorm, rankFIM, eigFIM] = LightSrcOptmLS(lightSrc, Phi0, targetL, targetPntFrame, targetPntNormals, reflectance, optOptions, intNoiseSigma)
% Fits a non-isotropic disk light source model of the given data in a
% non-linear least squares sense. Input sigma noise allows for the observability of
% the optimisation to be given.
% INPUTS:
%       lightSrc - point light source object
%       Phi0 - initial guess for unknown parameters 1x3 [max radiant intensity, mu, attenuation radius]
%       targetL - Measured irradiance
%       targetPntFrame - location of measured points on white reflectance target
%           w.r.t to the frame camera
%       targetPntNormals - normal at the measured location w.r.t to the
%           frame camera
%       reflectance - ratio of incoming light that is reflected in diffuse
%           manner for the current wavelength/frequency of light
%       optOptions - options struct for the lsqnonlin function
%       intNoiseSigma - zero-mean gaussian sigma noise in the measured irradiance
% OUTPUTS:
%       optPhi - optimised parameters 1x3 [max radiant intensity, mu, attenuation radius]
%       resNorm - sum of residuals after optimisation
%       rankFIM - rank of FIM
%       eigFIM - eigenvalues of FIM

%Author: Jasprabhjit Mehami, 13446277


%The optimisation function
f = @(H)LightSourceLeastSqrObj(H, targetL', targetPntFrame, ...
    targetPntNormals, lightSrc.get_ligSourLoc(), lightSrc.get_SourDirVec(), reflectance);

%Perform optimisation
[optPhi,resNormSqr] = lsqnonlin(f,Phi0, [], [], optOptions);

if any(optPhi < 0)
    optOptionsTemp = optimoptions('lsqnonlin', 'Algorithm', 'trust-region-reflective', 'SpecifyObjectiveGradient',true, 'CheckGradients', false, ...
        'MaxIterations', 1000000000, 'FunctionTolerance',1e-8, 'MaxFunctionEvaluations',1000000000, 'StepTolerance',1e-10, ...
        'FiniteDifferenceType', 'central');
    optOptionsTemp.Display = 'none';
    
    %Perform optimisation
    [optPhi,resNormSqr] = lsqnonlin(f,Phi0, [0,0,0], [], optOptionsTemp);
end




optPhi = optPhi';

if nargout < 2
    return;
end

resNorm = sqrt(resNormSqr);

if nargout < 3
    return;
end

%determine observability of NLS by calculating the Fisher Information
%Matrix
%Jacobian of radiance objective function w.r.t to unknown parameters
jac = JacobianLightRadiance(optPhi, targetPntFrame, ...
    targetPntNormals, lightSrc.get_ligSourLoc(), lightSrc.get_SourDirVec(), reflectance);

%create covariance matrix
cov = diag((1/(intNoiseSigma^2)).*ones(1,size(jac,1)));

FIM = (jac'/cov)*jac;

rankFIM = rank(FIM);
eigFIM = eig(FIM);
end

