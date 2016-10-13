% test empirical Bayes logistic regression code on simulated 1D example
addpath regressionTools

% 1.  Set up simulated example

% make 2D filter
nx = 16; % number of pixels
nt = 20; % number of frames 
nTrials = 2000; % number of trials

sig = 2;  % Template width
[xx,tt] = meshgrid(1:nx,1:nt);
wtsim = 2*exp(-((xx-nx/2).^2+(tt-1-nt/2).^2)/(2*sig^2)).*cos(xx-nx/2);
wts = wtsim(:)/5;
nw = nx*nt;
b = -1; % constant (DC term)

% Make stimuli & simulate response
stim = 1*(randn(nTrials,nw));
xproj = stim*wts+b;
pp = logistic(xproj);
yy = rand(nTrials,1)<pp;

% -- make plot ---
clf; 
subplot(221);
imagesc(wtsim);title('true filter (2D)');
iiw = 1:nw;
subplot(222);
plot(iiw,wts,'k');
title('true filter plotted as vector')
subplot(223);
xpl = min(xproj):.1:max(xproj);
plot(xproj,yy,'.',xpl,logistic(xpl), 'k');
xlabel('input'); ylabel('response');

errfun = @(w)(norm(w-wts).^2);  % error function handle

%% 2. Compute linear and mle estimates
xx = [stim, ones(nTrials,1)];  % regressors

% LS estimate
wls = xx\yy;  wlsplot = wls(1:nw)./norm(wls(1:nw))*norm(wts);


wmle=glmfit(xx(:,1:end-1), yy, 'binomial'); % remove bias parameter because glmfit has it's own
wmleplot=wmle(2:end); % remove bias parameter
wmle=[wmle(2:end); wmle(1)]; % reorder so the bias parameter is on the right. All the other code expects it that way
subplot(212);
plot(iiw,wts,'k',iiw,wlsplot,iiw,wmleplot);
legend('original', 'LS', 'glmfit');

%% 3. Compute Empirical Bayes logistic-regression estimate, AR1 2D prior

rhovals = 10.^(0:6)'; % grid over prior precision (hyperparameter)
avals = [.1 .5 .75 .9 .95 .99]'; % grid over correlation (AR1 hyperparameter)
rhoNull = .01;  % prior precision for other variables
[wRidge,hprsRidge,SDerrbarsRidge,HessRidge] = autoRegress_logisticRidge(xx,yy,[nt nx],rhoNull,rhovals,wmle);
[wAR1,hprsAR1,SDerrbarsAR1,HessAR1] = autoRegress_logisticAR1_2D(xx,yy,[nt nx],rhoNull,rhovals,avals,wmle);

plot(iiw,wts,'k'); hold on;
errorbar(iiw,wRidge(1:nw),2*SDerrbarsRidge(1:nw),'r');
errorbar(iiw,wAR1(1:nw),2*SDerrbarsAR1(1:nw),'c');
axis tight;
hold off;
legend('true', 'ridge', 'AR1');


% Err = [errfun(wAR1(1:nw))]

%% sparse prior in a smooth basis

addpath sparseglm/toolbox

nLevels=4;
step=.5;
fwhm=1.5; % full width of basis functions at half max


%Laplacian pyramid basis
B = get2DLaplacianPyramidBasis(nt,nx,nLevels,step,fwhm);
XB = stim*B; % project stimulus on the basis
U=ones(nTrials,1);
thefit = cvglmfitsparseprior(yy,XB,U,getcvfolds(length(yy),5));
imestlp2 = reshape(B*thefit.w,nt,nx);

%% plot the different methods
figure(121); clf
subplot(161)
imagesc(wtsim); title('true filter')
subplot(162)
imagesc(reshape(wlsplot, [nt nx])); title('linear')
subplot(163)
imagesc(reshape(wmleplot, [nt nx])); title('glmfit')
subplot(164)
imagesc(reshape(wRidge(1:nw), [nt nx])); title('ridge')
subplot(165)
imagesc(reshape(wAR1(1:nw), [nt nx])); title('AR1')
subplot(166)
imagesc(imestlp2); title('Sparse and smooth')


Errs = [errfun(wRidge(1:nw)), errfun(wAR1(1:nw)) errfun(imestlp2(:))]
