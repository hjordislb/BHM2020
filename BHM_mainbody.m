%***********************************************************************
%  A Spatial Bayesian Hierarchical model (BHM) for ground motion amplitudes 
% 
%  PURPOSE : A new Bayesian hierarchical model (BHM) for the variability of
%            PGA is developed across a small-aperture array. The proposed 
%            model incorporates many of the commonly used seismic parameters
%            and offers a flexible probabilistic framework for multilevel 
%            modeling of PGA that accounts for the effects of the earthquake
%            source, propagation path effects and localized site effects, 
%            along with their respective variabilities. 
%            MCMC algorithm is used to sample from the posterior density of 
%            the proposed BHM.
%  
%  INPUT : - Ground-motion amplitudes 
%          - inter-station distance or station locations 
%          - seismic parameters (local magnitude, source-to-site distance, 
%                                depth, direction)
%  
%  OUTPUT : Posterior Distribution of hyperparameters and latent parameters
%           along with statistical convergence diagnostics
% 
%  WRITTEN : Sahar Rahpeyma, Ph.D.      
%            Faculty of Civil & Environmental Engineering 
%            Earthquake Engineering Research Centre                   
%            University of Iceland                                    
%            Email: sahar@hi.is                                       
%            Ref: https://onlinelibrary.wiley.com/doi/abs/10.1002/env.2497
%                 Rahpeyma et al. (2018) Bayesian hierarchical model for 
%                 variations in earthquake peak ground acceleration within 
%                 small-aperture arrays, Environmetrics, 29(3), e2497 
%                                                                       
%  DATE    : May.2020                    
%***********************************************************************
%% ================
% INPUT: 
% =================
clear;
FS = filesep;
% NOTE : here you need to define the path to your file locations as below:
MAIN = 'C:\Users\hjord\Documents\Sumarvinna-2020-HI\ProjectCode\Sahar_Rahpeyma_BHM2020'; %SR 
DAT  = '\Data';

% If input data is saved in .m files use:
load ([MAIN FS DAT FS 'PGA_Obs']);  % Ground motion amplitudes (Mat)
load ([MAIN FS DAT FS 'dij']);      % Inter-station distance (Mat)
load ([MAIN FS DAT FS 'BAz']);      % Back-azimuth i.e. direction (Mat)
load ([MAIN FS DAT FS 'M']);        % Local magnitude (Vec)
load ([MAIN FS DAT FS 'Depth']);    % Depth (Vec)
load ([MAIN FS DAT FS 'R_HYP']);    % Hypocentral distance (Mat)
load ([MAIN FS DAT FS 'R_EPI']);    % Epicentral distance (Mat)


% If your input data is saved in .txt, .dat, or .csv, .xls, .xlsb, .xlsm,
% .xlsx, .xltm, .xltx, or .ods files use:

% PGA_Obs = readmatrix([MAIN FS DAT FS 'PGA_Obs.xlsx']);  % Ground motion amplitudes (M)
% dij = readmatrix ([MAIN FS DAT FS 'dij.xlsx']);      % Inter-station distance (M)
% BAz = readmatrix ([MAIN FS DAT FS 'BAz.xlsx']);      % Back-azimuth i.e. direction (M)
% M = readmatrix ([MAIN FS DAT FS 'M.xlsx']);        % Local magnitude (V)
% Depth = readmatrix ([MAIN FS DAT FS 'Depth.xlsx']);    % Depth (V)
% R_HYP = readmatrix ([MAIN FS DAT FS 'R_HYP.xlsx']);    % Hypocentral distance (M)
% R_EPI = readmatrix ([MAIN FS DAT FS 'R_EPI.xlsx']);    % Epicentral distance (M)

%% ==================================
% DEFINE SOME GRAPHICS PARAMETERS: 
% ===================================

FigPath = [MAIN FS 'Figs'];
figtype = '-djpeg'; 
figdpi  = '-r400'; 
LS = 'LineStyle';
fs = 'fontsize';
CL = 'color'; 
PO = 'position';
FW = 'FontWeight';
ROT = {'rotation',90; 'rotation',0};
MFC = 'markerfacecolor';
MS  = {'markersize',5; 'markersize',10;
      'markersize',15;'markersize',20};
MEC = 'MarkerEdgeColor';                 
LW  = {'LineWidth',1.2;'LineWidth',1.5;  
      'LineWidth',2.0;'LineWidth',3.0};
LST = {'-'; '--'; ':'};                   
clr = {'k'; 'r'; 'b'; 'm'; 'c'};          
UnitNorm = {'Units','Normalized'};         
%% ================================
%PREPARE INPUT FOR MODELING BHM:
% =================================
y=PGA_Obs;

y (y==0)=-9999; % missed data (i.e. stations without records) denoted by -9999
NS = size(y,2); % sites: 10   (deltaS2S)
NT = size(y,1); % event: 60   (deltaB) note: selected from original dataset
NY = numel(y);  % all data: [10 * 60] = 600 (without missing points)
Y = reshape(y',NY,1);        % vector: [all_sites(time 1); time 2...]
Dis = reshape(R_HYP',NY,1);  % Vector of Distance velues
BA  = reshape(BAz',NY,1);    % Vector of Back-Azimuth values 
for i=1:NY
    if   Y(i)~= -9999
         Y(i) = log10(Y(i)); 
    else
         Y(i) = -9999; 
    end
end
Z1 = kron(eye(NT),ones(NS,1));  % (A.2)
Z2 = kron(ones(NT,1),eye(NS));  % (A.3)
X = [ones(NT,1) M log10(R_HYP(:,1)) Depth]; 
Z1(Y==-9999,:) = 0;
Z2(Y==-9999,:) = 0;
Y(Y==-9999) = 0;
ZX = Z1*X;
DEP = Z1*Depth;
zp = find(Y~=0);   % find and exclude missing data 
NNZY = length(zp); % total number of observations 
Y = Y(zp);
Z1 = Z1(zp,:);   
Z2 = Z2(zp,:);
ZX = ZX(zp,:);
Dis = Dis(zp);
DEP = DEP(zp);
BA = BA(zp);
ZX(:,3) = log10(Dis);
ZX(:,4) = DEP;   % ZX matrix is the same as X matrix 
W = [ZX Z1 Z2];  % Eq.(A.4): W is the same as K matrix in Rahpeyma et al. (2018)
%% =========================================== %
% Hyperparameters and latent parameters Prior:
% ============================================ % 
Lam_phi = ceil(sqrt(var(Y)/4)^-1); 
Delta = 0.7; 
lam_theta = [Lam_phi  Lam_phi  Lam_phi  Delta  Lam_phi  Delta];
% Beta    = [B1   B2   B3  B4 ] 
mu_beta   = [0    0    0    0  ]';
var_beta  = [100  100  100  100].^2';
% parameters and prior pdf
parm_hyp = {'\phi_R',...
            '\tau',...
            '\phi_{SS}','\Delta_{SS}',...
            '\phi_{S2S}','\Delta_{S2S}'};
prior_hyp = repmat({'exppdf'},1,length(lam_theta));
parm_lat = {'beta','delta B','delta S2S'};
prior_lat= repmat({'normpdf'},1,3);


  

%% ===================================== %% 
% Global mode optimization: 
% ====================================== %

msx_a = [ 0.01  0.01  0.01  0.100  0.01  0.01 ];
msx_b = [ 0.99  0.99  0.99  1.00   0.99  0.99 ];
Ngen = 1000; %  
Nbest = 50;  % closest values to the mode 
Niter = 1;  % itterations of uniform dist.
Npar =numel(msx_a); 
% Generating from uniform distirbution 
thetax = repmat(msx_a,Ngen,1) + ...
repmat(msx_b-msx_a,Ngen,1).*rand(Ngen,Npar); 
thetax (:,6)= 0.06;
% Generating from normal distribution     
for b=1:Ngen
    PU(b) = BHM_lnpost(...
           thetax(b,:),lam_theta,Y,W,dij,mu_beta,var_beta,zp);
end

Ir=PU==real(PU);
pR=PU(Ir);
[~,Is]=sort(pR);
best50=Is(end-49:end); % find best (i.e. highest) loglik.values
msxR  = thetax(Ir,:);
PU50  = pR(best50);
msx50 = msxR(best50,:);
%mean and standard deviation of the best 50 samples 
msx_m = mean(msx50);
msx_s = std(msx50);

eigen_vals = [-1,-1,-1,-1,-1];
fprintf('Finding postitive eigenvalues... \n')
while(any(eigen_vals<0))
    for i=1:Niter
        % Generating normal dist. from uniform distirbution 
        msx = repmat(msx_m,Ngen,1) + repmat(msx_s,Ngen,1).*abs(randn(Ngen,Npar));
        msx(:,6) = 0.06;
        PN = BHM_lnpost(msx ,lam_theta,Y,W,dij,mu_beta,var_beta,zp);
        [~,maxi] = max(PN);
        mode_msx = msx(maxi,:); % MAXIMUM OF P obtained by these three param.
    end
    mode_theta = mode_msx; % it has to be checked for positive COV
    
                     
    h = [1 1 1 1 1 1]*0.0001;     % step size numerical Hessian i.e. 

    H = my_hessian(@BHM_lnpost,mode_theta, h,...       % Generalised hessian Matrix Eq.(6)
                    lam_theta,Y,W,dij,mu_beta,var_beta,zp); 
    H = H(1:5,1:5);                  % delta S2S is exclued            
    NPhyp = length(parm_hyp);        % number of hyperparameters 
    scale = 2.38^2/(NPhyp-1);        % scale parameter 
    Pcov = -scale*(H\eye(NPhyp-1));  % precision matrix (Roberts et al.,1997)
    eigen_vals = eig(Pcov)           % checking if Pcov is positive definite
end

%% =======================================
% Fixing the mode:
% ========================================
% When a suitable mode has been found, it
% can be fixed here as "mode_theta"

mode_theta = mode_theta
h = [1 1 1 1 1 1]*0.0001;     % step size numerical Hessian
H = my_hessian(@BHM_lnpost,mode_theta, h,...       % Hessian Matrix Eq.(6)
                    lam_theta,Y,W,dij,mu_beta,var_beta,zp); 
H = H(1:5,1:5);                  % delta S2S is exclued            
NPhyp = length(parm_hyp);        % number of hyperparameters 
scale = 2.38^2/(NPhyp-1);        % scale parameter 
Pcov = -scale*(H\eye(NPhyp-1));  % precision matrix (Roberts et al.,1997)
eigen_vals = eig(Pcov)           % checking if Pcov is positive definite


%% =======================================
% Gibbs sampler: Setup
% ========================================
% memory initialization:

NT = 1000;      % total number of iterations (For accurate results, try 10.000-50.000)
NL = 0.75*NT;    % sample size L to be drawn (after burn-in)
NB = 0.25*NT;    % burn-in samples at beginning 
NC = 4;          % number of chains to compute
NT = NL + NB;  
Nt = NY/NS;      % number of events
NPhyp = length(parm_hyp);    % number of hyperparameters 
NBet  = length(mu_beta);     % Beta coefficients (\Beta)
NPlat = NBet+NS+Nt;          % number of latent parameters 
% prepare memory for parameters:
% Phyp = nan(NPhyp,NC,NT); 
Phyp = cell(NC,1);          % hyperparameters 
Phyp(:) = {nan(NT,NPhyp)};  
Plat = cell(NC,1);          % latent parameters
Plat(:) = {nan(NT,NPlat)};  
LPhyp = cell(NC,1);         % hyperparameetrs [log-scale]
LPhyp(:) = {nan(NT,NPhyp)}; 
% additional constants:
%scale = 2.38^2/(NPhyp-1);          
%Pcov  = scale*(-H\eye(NPhyp-1));   % covariance of proposal pdf (Roberts et al.,1997)
Pchol = chol(Pcov);
MU_eta = [mu_beta;zeros(NS+Nt,1)]; % latent prior mean
cov_beta = diag(var_beta);         % latent prior covariance
IL = eye(NPlat);
IY = eye(NNZY);

%% =======================================
% Gibbs sampler loop  (section 3.2 | Posterior inference)
% =======================================
% Prepare variables to SAVE the posterior Model Parameters in 
% sim_hyp: hyper parameters & sim_latent: Latent parameters
sim_hyp = sprintf('%s\\mat\\sim_hyper.mat',MAIN);
sim_lat = sprintf('%s\\mat\\sim_latent.mat',MAIN);

% How many Loops? (Optional)
% matlabpool open 4 : Matlab Versions 2007-2014
parpool(2)        % : Matlab Versions after 2014
for c = 1:NC   
    
    % INITIAL values for all chains and parameters (hyper & latent)
    % -------------------- Hyperparameters -------------------- % 
    Phyp{c}(1,1:5) = mvnrnd(mode_theta(1:5), Pcov); % generate the first sample using mode values (in log)
    Phyp{c}(1,6) = 0.06;                       % fix the range parameter of site effect
%     Phyp{c}(1,:) = exp(LPhyp{c}(1,:));               % hyperparameters   
    % log posterior p(theta|y)
    p1 = BHM_lnpost(Phyp{c}(1,:),lam_theta,Y,W,dij,mu_beta,var_beta,zp);          
    % -------------------- Latent parameters ------------------ % 
    cov_S2S = exp(2*log(Phyp{c}(1,5)))*exp(-dij/exp(log(Phyp{c}(1,6))));   % Eq.(A.8) 
    
    cov_SS_sub = exp(2*log(Phyp{c}(1,3))).*(1+dij/exp(log(Phyp{c}(1,4)))).*...
                 exp(-dij/exp(log(Phyp{c}(1,4))));  %Eq.(A.6)
                              
    cov_SS = kron(eye(Nt,Nt),cov_SS_sub); 
    cov_SS = cov_SS(zp,zp);
    cov_tau = exp(2*log(Phyp{c}(1,2)))*eye(Nt);
   
    % latent prior precision  Eq.(A.16)
    cov_eta = [cov_beta        zeros(NBet,Nt)   zeros(NBet,NS);   
               zeros(Nt,NBet)  cov_tau          zeros(Nt,NS);
               zeros(NS,NBet)  zeros(NS,Nt)     cov_S2S];
    SEI = cov_eta\IL; 
    % Cov_eta inv.
    SYI = (exp(2*log(Phyp{c}(1,1)))*IY + cov_SS)\IY; % data error precision (Sig_y inv.)     Eq.(A.18)
    SLP = (SEI + W'*SYI*W)\IL;                       % Sigma{lat,post} latent posterior cov. Eq.(A.19)
    Plat{c}(1,:) = (SLP*(SEI*MU_eta + W'*SYI*Y) + ...
                   (chol(SLP))'*randn(NPlat,1))';    % Eq.(A.20)
               
    % ***************************************************************** %
    % PROPOSE NEW VALUES for all chains and parameters (hyper & latent)
    % ***************************************************************** %
    for a = 2:NT 
    % -------------------- Hyperparameters -------------------- %  
    Phyp{c}(a,1:5) = mvnrnd(Phyp{c}(a-1,1:5), Pcov);
    Phyp{c}(a,6)   = 0.06;

%         % get probability ...
          p2 = BHM_lnpost(Phyp{c}(a,:),lam_theta,Y,W,dij,...
                mu_beta,var_beta,zp);
%         % get probability ...
%         p2 = BHM_lnpost(Phyp{c}(a,:),lam_theta,Y,W,dij,...
%                 mu_beta,var_beta,zp); 
%         % acceptance or rejection:  Eq.(8) & Eq.(9)
        % REJECT if ...  
        if rand > exp(p2 - p1) 
            Phyp{c}(a,:) = Phyp{c}(a-1,:);
%             LPhyp{c}(a,:) = LPhyp{c}(a-1,:);
            Plat{c}(a,:) = Plat{c}(a-1,:);
            
        else 
    % -------------------- Latent parameters ------------------ %         
        cov_S2S = exp(2*log(Phyp{c}(a,5)))*exp(-dij/exp(log(Phyp{c}(a,6))));

        cov_SS_sub = exp(2*log(Phyp{c}(a,3))).*(1+dij/exp(log(Phyp{c}(a,4)))).*...
                exp(-dij/exp(log(Phyp{c}(a,4))));

        cov_SS = kron(eye(Nt,Nt),cov_SS_sub); 
        cov_SS = cov_SS(zp,zp);
        cov_tau = exp(2*log(Phyp{c}(a,2)))*eye(Nt);

        % latent prior precision
        cov_eta = [cov_beta        zeros(NBet,Nt)   zeros(NBet,NS);   
                   zeros(Nt,NBet)  cov_tau          zeros(Nt,NS);
                   zeros(NS,NBet)  zeros(NS,Nt)     cov_S2S ];
        SEI = cov_eta\IL; % Cov_eta inv.      
        SYI = (exp(2*log(Phyp{c}(a,1)))*IY + cov_SS)\IY;     % data error precision (Sig_y inv.)
        SLP = (SEI + W'*SYI*W)\IL;                           % Sigma{lat,post} latent posterior cov.
        Plat{c}(a,:) = (SLP*(SEI*MU_eta + W'*SYI*Y) + ...
                       (chol(SLP))'*randn(NPlat,1))';
        % Replace new probability
        p1 = p2;
        end
        fprintf('iter. %d/%d, chain %d/%d done\n',a,NT,c,NC); % display (iteration : xx (%Chain xx))
    end
end
% matlabpool close
delete(gcp)  

% save simulated latent and hyperparameters:
save(sim_hyp,'Phyp');
save(sim_lat,'Plat');

%% ===================================================== % 
% Reshape simulated parameters obtained by Gibbs sampler:
% ====================================================== % 
NN = NL*NC;
parh  = nan(NT,NC,NPhyp);
parhc = nan(NL*NC,NPhyp);
parl  = nan(NT,NC,NPlat);
parlc = nan(NL*NC,NPlat);

for b = 1:NPhyp
    for a = 1:NC
        parh(:,a,b) = Phyp{a}(:,b);
    end
    parhc(:,b) = reshape(parh(NB+1:NT,:,b),NN,1);
end

for b = 1:NPlat
    for a = 1:NC
        parl(:,a,b) = Plat{a}(:,b);
    end
    parlc(:,b) = reshape(parl(NB+1:NT,:,b),NN,1);
end

% statistics: [mean, std, 2.5%, 16%, 25%, 50%, 75%, 84%, 97.5%]
sort_hyp = sort(parhc,1);
sort_lat = sort(parlc,1);
ii = round([0.025 0.16 0.25 0.5 0.75 0.84 0.975]*NN);
% ***************************************************************** % 
% Mean, SD, and Percentile for hyperparameters and latent parameters  
% ***************************************************************** %
%          [mean SD 2.5% 16% 25% 50% 75% 84% 97.5%]' x NPhyp hyperparameters
hyp_stat = [mean(sort_hyp,1); std(sort_hyp,0,1); sort_hyp(ii,:)];
%          [mean SD 2.5% 16% 25% 50% 75% 84% 97.5%]' x (NBet+NT+NS) latent parameters
lat_stat = [mean(sort_lat,1); std(sort_lat,0,1); sort_lat(ii,:)];

%% ===================================================== % 
% CONVERGENCE DIAGNOSTICS : 
% TRACEPLOTS, HISTOGRAMS, GELMAN-RUBIN, AUTOCORRELATION STATISTICS
% ===================================================== % 
% ********************************** % 
% Gelman-Rubin and Auto-correlation 
% ********************************** % 
II = 50:100:NT; %could change steps from 10 to 100 (or another value)
ii = numel(II); 
R = nan(NPhyp+NPlat,ii); Neff = R; lag1 = R; acptR = R;

% hyperparameters (\theta)
for b = 1:NPhyp-1
    for i = 1:ii
        k = II(i);
        [R(b,i),Neff(b,i),lag1(b,i),acptR(b,i)] = gpar(parh(1:k,:,b));
    end
end
% latent parameters : constant coefficients (\beta)  
for a=1:NBet 
    for i=1:ii;
        k = II(i);
        [R(a+NPhyp,i),Neff(a+NPhyp,i),lag1(a+NPhyp,i),acptR(a+NPhyp,i)] =...
            gpar(parl(1:k,:,a));
    end
end
% latent parameters : event terms (\delta B)  
for a=1:Nt 
    for i=1:ii;
        k = II(i);
        [R(a+NPhyp+NBet,i), Neff(a+NPhyp+NBet,i), lag1(a+NPhyp+NBet,i), ...
            acptR(a+NPhyp+NBet,i)] = gpar(parl(1:k,:,a));
    end
end
% latent parameters: site-terms (\delta-S2S)
for a=1:NS 
    for i=1:ii;
        k = II(i);
        [R(a+NPhyp+Nt+NBet,i),Neff(a+NPhyp+Nt+NBet,i),...
            lag1(a+NPhyp+Nt+NBet,i),acptR(a+NPhyp+Nt+NBet,i)] =...
            gpar(parl(1:k,:,a));
    end
end
% ********************************** % 
% lag-autocorrelation 
% ********************************** % 
mlag = 50;
aut_lag = zeros(NPhyp+NPlat,mlag+1);
for a=1:NPhyp  % \theta : Hyperparameters
    for b=1:NC
        XCR = xcorr(parh(NB+1:end,b,a)-mean(parh(NB+1:end,b,a)), ...
            mlag,'coeff');
        aut_lag(a,:) = aut_lag(a,:) + XCR(mlag+1:end)'/NC;
    end
end
for a=1:NBet % \beta : constant coefficients of log(mu_{es}) 
    for b=1:NC
        XCR = xcorr(parl(NB+1:end,b,a) - mean(parl(NB+1:end,b,a)), ...
            mlag,'coeff');
        aut_lag(a+NPhyp,:) = aut_lag(a+NPhyp,:) + XCR(mlag+1:end)'/NC;
    end
end
for a=1:NS % \DeltaS2S : site-terms 
    for b=1:NC
        XCR = xcorr(parl(NB+1:end,b,a)-mean(parl(NB+1:end,b,a)), ...
            mlag,'coeff');
        aut_lag(a+NPhyp+Nt+NBet,:) = aut_lag(a+NPhyp+Nt+NBet,:) + ...
            XCR(mlag+1:end)'/NC;
    end
end
% ********************************** % 
% TRACEPLOTS 
% ********************************** % 
%% HYPERPARAMETERS
close all
SC = get(0,'ScreenSize');
fig1 = figure;

yl=[0 0.10 ;  %\phi_R
    0 0.25 ;  %\tau
    0 0.30 ;  %\phi_SS
    0 2.00 ;  %\Delta_SS
    0 0.3];  %\phi_S2S
h  = 0.15;
y0 = 0.83;
dy = 0.18;
for i = 1:NPhyp-1
    % ------------------------------------- % 
    % traceplot for all the chains 
    % ------------------------------------- % 
    sp = subplot('position', [0.1 y0-(i-1)*dy 0.45 h]);
    plot(squeeze(parh(:,:,i))); 
    hold on
    plot(squeeze(parh(1:NB,:,i)),'color', [0.5 0.5 0.5]);
    Ylab = ylabel(sprintf('%s',parm_hyp{i}), fs, 12);
    set (Ylab, ROT{2,:}, UnitNorm{1,:}, PO, [-.15 0.40 0]);
    if i~=NPhyp-1
        set(gca, 'xticklabel', [])
    else
        xlabel('Total iteration', 'fontsize', 12);
    end
        ylim(yl(i,:));
   % ------------------------------------- %  
   % histograms, rotated 90 [deg]
   % ------------------------------------- %  
    subplot('position', [0.56 y0-(i-1)*dy 0.105 h]);
    [ybin,xval] = hist(sort_hyp(:,i), 50);
    b = barh(xval, (ybin./max(ybin)), 'FaceColor', [0.20 0.5 0.9],...
        'EdgeColor', [0.5 0.5 0.5]); hold on   
    set(gca,'yticklabel',[]); 
    set(gca,'xticklabel',[]);
    box off; axis off
    ax=axis;
    plot(ax(1:2), hyp_stat(1,i)*[1 1]', 'r-');
    plot(ax(1:2), hyp_stat(4,i)*[1 1]', 'r--'); 
    plot(ax(1:2), hyp_stat(8,i)*[1 1]', 'r--');
    ylim(yl(i,:));
   % ------------------------------------- %  
   % Gelman-Rubin
   % ------------------------------------- %  
   subplot('position', [0.68 y0-(i-1)*dy .125 h]);
   plot(II,R(i,:))
    if i~=NPhyp-1
        set(gca,'xticklabel', [])
    else
        xlabel('Total iteration', 'fontsize',12);
    end
    ylim([1 2 ])
    xlim([0 NT])
    hline = refline([0 1.1]); 
    set(hline, CL, 'k', 'LineStyle', LST{3}, LW{2,:})
    if i==1
        TX = text(0.1, 0.1 ,'Gelman-Rubin',fs,10);
        set(TX, UnitNorm{1,:}, PO, [0.05,0.88,0]);
    end
   % ------------------------------------- %  
   % Autocorrelation
   % ------------------------------------- % 
   subplot('position', [0.85 y0-(i-1)*dy 0.125 h]); 
   bar(1:mlag, aut_lag(i,1:end-1), 'FaceColor',[0.5 0.5 0.5],...
       'EdgeColor',[0.5 0.5 0.5]);
   if i~=NPhyp-1
       set(gca, 'xticklabel', []);
   else
        xlabel('Lag', fs, 12); 
   end
   if i==1
        TX = text(0.1, 0.1 ,'AutoCorr', fs, 10);
        set(TX, UnitNorm{1,:}, PO, [0.05,0.88, 0]);
    end
   xlim([0 50])
   ylim([-0.00 1])
end
% save figure
spr = sprintf('%s\\traceplot_hpar.jpg',FigPath); 
print(figtype , figdpi , spr);

%% LATENT PARAMETERS
% \Beta
fig2 = figure;
yl=[0 1.5;0.5 1;-5 -2.0; 0 0.2];
y0 = 0.80; 
dy = 0.23; 
h  = 0.185;
for i = 1:NBet
    % ------------------------------------- % 
    % traceplot for all the chains 
    % ------------------------------------- % 
    sp = subplot('position',[0.1 y0-(i-1)*dy 0.45 h]);
    plot(squeeze(parl(:,:,i))); hold on
    plot(squeeze(parl(1:NB,:,i)),'color',[0.5 0.5 0.5]);
    Ylab = ylabel(sprintf('\\%s_%d', parm_lat{1},i),fs,12);
    ylim(yl(i,:));
    set (Ylab, ROT{2,:}, UnitNorm{1,:}, PO, [-0.13 0.35 0]);
    if i~=NBet
        set(gca,'xticklabel',[])
    else
        xlabel('Total iteration', fs,12);
    end
   % ------------------------------------- %  
   % histograms, rotated 90 [deg]
   % ------------------------------------- %  
    subplot('position',[0.56 y0-(i-1)*dy .105 h]);
    [ybin,xval] = hist(sort_lat(:,i),70);
    barh(xval,(ybin/max(ybin)),'FaceColor',[0.20 0.5 0.9],...
        'EdgeColor',[0.5 0.5 0.5]); hold on 
    ylim(yl(i,:))
    ax = axis;
    plot(ax(1:2), lat_stat(1,i)*[1 1]', 'r-' ); % mean posterior
    plot(ax(1:2), lat_stat(4,i)*[1 1]', 'r--'); % 16% posterior percentile 
    plot(ax(1:2), lat_stat(8,i)*[1 1]', 'r--'); % 84% posterior percentile
    set(gca, 'yticklabel', []); 
    set(gca, 'xticklabel', []);
    box off; axis off
   % ------------------------------------- %  
   % Gelman-Rubin
   % ------------------------------------- %  
   subplot('position',[0.68 y0-(i-1)*dy 0.125 h]);
   plot(II,R(NPhyp+i,:))
    if i~=NBet
        set(gca,'xticklabel',[])
    else
        xlabel('Total iteration', fs,12);
    end
    ylim([1 2 ])
    xlim([0 NT])
    hline = refline([0 1.1]); 
    set(hline,'Color','k','LineStyle',':','linewidth',1.5)
    if i==1
        TX =text(0.1, 0.1 ,'Gelman-Rubin',fs,10);
        set(TX,UnitNorm{1,:},'Position', [0.05,0.88, 0.0]);
    end
   % ------------------------------------- %  
   % Autocorrelation
   % ------------------------------------- % 
   subplot('position',[0.85 y0-(i-1)*dy .125 h]); 
   bar(1:mlag, aut_lag(NPhyp+i,1:end-1), 'FaceColor', [0.5 0.5 0.5],...
       'EdgeColor', [0.5 0.5 0.5]);
   if i~=NBet
       set(gca, 'xticklabel', [])
   else
        xlabel('Lag', fs, 12), 
   end
   if i==1
        TX =text(0.1, 0.1 , 'AutoCorr', fs, 10);
        set(TX,UnitNorm{1,:},'Position', [0.05,0.88, 0]);
    end
    ylim([-0.05 1])
    xlim([0.00 50])
end

% save figure
spr = sprintf('%s\\traceplot_Beta.jpg',FigPath); 
print(figtype , figdpi , spr);

%% \delta_S2S
fig3 =  figure;
y0 = 0.920; 
dy = 0.091; 
h  = 0.072;
Scode = {'IS601'; 'IS602'; 'IS603'; 'IS604'; 'IS605'; 'IS607'; 'IS608'; ...
         'IS609'; 'IS611'; 'IS612'};
for s=1:NS
    % ------------------------------------- %  
    % traceplot for all the chains 
    % ------------------------------------- %  
    sp = subplot('position',[0.1 y0-(s-1)*dy 0.45 h]);
    plot(squeeze(parl(:,:,Nt+NBet+s))); hold on
    plot(squeeze(parl(1:NB,:,Nt+NBet+s)),'color',[0.5 0.5 0.5]);
    Ylab = ylabel(sprintf('\\%s_{%d}', parm_lat{3},s), fs, 12);
    ylim([-0.25 0.25])
    set(gca,'YTick',[-0.25  0.25]);
    set(gca,'YTickLabel',[-0.25  0.25])
    set (Ylab, 'rotation',0, UnitNorm{1,:}, PO, [-0.13 0.40 0.00]);
    if s~=NS
        set(gca,'xticklabel',[])
    else
        xlabel('Total iteration', fs,12);
    end
    TX =text(0.1, 0.1 , sprintf('%s', Scode{s}));
    set(TX, UnitNorm{1,:}, PO, [0.02, 0.88, 0.00]);
    Zline = refline([0 0]); 
    set(Zline, CL, clr{1}, LS, LST{2}, LW{1,:})
   % ------------------------------------- %  
   % histograms, rotated 90 [deg]
   % ------------------------------------- %  
    subplot('position',[0.56 y0-(s-1)*dy 0.105 h]);
    [ybin,xval] = hist(sort_lat(:,Nt+NBet+s),50);
    b = barh(xval,(ybin/max(ybin)),'FaceColor',[0.20 0.5 0.9],...
        'EdgeColor',[0.5 0.5 0.5]); hold on  
    ylim([-0.25 0.25])
    ax=axis;
    plot(ax(1:2),lat_stat(1,Nt+NBet+s)*[1 1]','r-');  % mean posterior
    plot(ax(1:2),lat_stat(4,Nt+NBet+s)*[1 1]','r--'); % 16% posterior percentile
    plot(ax(1:2),lat_stat(8,Nt+NBet+s)*[1 1]','r--'); % 84% posterior percentile 
    set(gca,'yticklabel',[]); 
    set(gca,'xticklabel',[]);
    box off; axis off
   % ------------------------------------- %  
   % Gelman-Rubin
   % ------------------------------------- %  
   subplot('position',[0.68 y0-(s-1)*dy 0.125 h]);
   plot(II,R(NPhyp+Nt+NBet+s,:))
    if s~=NS
        set(gca,'xticklabel',[])
    else
        xlabel('Total iteration', fs,12);
    end
    ylim([1 2])
    xlim([0 NT])
    hline = refline([0 1.1]); 
    set(hline,'Color','k','LineStyle',':','linewidth',1.5);
    if s==1
        TX = text(0.1, 0.1 , 'Gelman-Rubin', fs, 10);
        set(TX,UnitNorm{1,:}, 'Position', [0.05,0.88, 0]);
    end
   % ------------------------------------- %  
   % Autocorrelation
   % ------------------------------------- % 
   subplot('position',[0.85 y0-(s-1)*dy .125 h]); 
   bar(1:mlag, aut_lag(NPhyp+Nt+NBet+s,1:end-1),'FaceColor',[0.5 0.5 0.5],...
       'EdgeColor', [0.5 0.5 0.5]);
   if s~=NS
       set(gca, 'xticklabel',[])
   else
        xlabel('Lag', fs, 12), 
   end
   if s==1
        TX =text(0.1, 0.1 ,'AutoCorr', fs, 10);
        set(TX, UnitNorm{1,:}, PO, [0.05,0.88, 0]);
    end
    ylim([-0.05 1 ])
    xlim([ 0.00 50])
    % 
end
% save figure
spr = sprintf('%s\\traceplot_deltaS2S.jpg',FigPath); 
print(figtype , figdpi , spr);

