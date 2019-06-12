% replicates matlab's hmc example to test for understanding

% true parameters (I would know this for simulated data (which I should
% test this on before experimental data)
NumPredictors = 2;
trueIntercept = 2;
trueBeta = [3;0];
trueNoiseSigma = 1;

% data that fit will be performed on
NumData = 100;
X = rand(NumData,NumPredictors);
mu = X*trueBeta + trueIntercept;
y = normrnd(mu,trueNoiseSigma);

% set up priors
InterceptPriorMean = 0;
InterceptPriorSigma = 10;
BetaPriorMean = 0;
BetaPriorSigma = 10;
LogNoiseVarianceMean = 0;
LogNoiseVarianceSigma = 2;

% define logpdf function
logpdf = @(Parameters)logPosterior(Parameters,X,y, ...
    InterceptPriorMean,InterceptPriorSigma, ...
    BetaPriorMean,BetaPriorSigma, ...
    LogNoiseVarianceMean,LogNoiseVarianceSigma);

% create HMC Sampler
Intercept = randn;
Beta = randn(NumPredictors,1);
LogNoiseVariance = randn;
startpoint = [Intercept;Beta;LogNoiseVariance];
smp = hmcSampler(logpdf,startpoint,'NumSteps',50);

% estimates maximum a-posteriori 
[MAPpars,fitInfo] = estimateMAP(smp,'VerbosityLevel',1);
MAPIntercept = MAPpars(1);
MAPBeta = MAPpars(2:end-1);
MAPLogNoiseVariance = MAPpars(end);

% tunes sample parameters
[smp,tuneinfo] = tuneSampler(smp,'Start',MAPpars);

% log posterior probability for the function
function [logpdf, gradlogpdf] = logPosterior(Parameters,X,Y, ...
    InterceptPriorMean,InterceptPriorSigma, ...
    BetaPriorMean,BetaPriorSigma, ...
    LogNoiseVarianceMean,LogNoiseVarianceSigma)


% Unpack the parameter vector
Intercept        = Parameters(1);
Beta             = Parameters(2:end-1);
LogNoiseVariance = Parameters(end);
% Compute the log likelihood and its gradient
Sigma                   = sqrt(exp(LogNoiseVariance));
Mu                      = X*Beta + Intercept;
Z                       = (Y - Mu)/Sigma;
loglik                  = sum(-log(Sigma) - .5*log(2*pi) - .5*Z.^2);
gradIntercept1          = sum(Z/Sigma);
gradBeta1               = X'*Z/Sigma;
gradLogNoiseVariance1	= sum(-.5 + .5*(Z.^2));
% Compute log priors and gradients
[LPIntercept, gradIntercept2]           = normalPrior(Intercept,InterceptPriorMean,InterceptPriorSigma);
[LPBeta, gradBeta2]                     = normalPrior(Beta,BetaPriorMean,BetaPriorSigma);
[LPLogNoiseVar, gradLogNoiseVariance2]  = normalPrior(LogNoiseVariance,LogNoiseVarianceMean,LogNoiseVarianceSigma);
logprior                                = LPIntercept + LPBeta + LPLogNoiseVar;
% Return the log posterior and its gradient
logpdf               = loglik + logprior;
gradIntercept        = gradIntercept1 + gradIntercept2;
gradBeta             = gradBeta1 + gradBeta2;
gradLogNoiseVariance = gradLogNoiseVariance1 + gradLogNoiseVariance2;
gradlogpdf           = [gradIntercept;gradBeta;gradLogNoiseVariance];
end

% log prior probability function
 function [logpdf,gradlogpdf] = normalPrior(P,Mu,Sigma)
Z          = (P - Mu)./Sigma;
logpdf     = sum(-log(Sigma) - .5*log(2*pi) - .5*(Z.^2));
gradlogpdf = Z./Sigma;
end
