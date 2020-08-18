function p = BHM_lnpost_2nd(theta,lam_theta,Y,W,dij,mu_beta,var_beta,zp)

% This Routine calculate log-posterior of hyper parameters p(theta|y) 
% for Bayesian Hierarchical Model (BHM) of Ground motion Amplitude 

%  WRITTEN : Sahar Rahpeyma, Ph.D.      
%            Faculty of Civil & Environmental Engineering 
%            Earthquake Engineering Research Centre                   
%            University of Iceland                                    
%            Email: sahar@hi.is                                       
%            Ref: Rahpeyma et al. (2018) 
%                                                                       
%  DATE    : May.2020               VERSION:  3.0             

% -------------------------------------------------------------------------
% INPUT: Y = data column vector 
%        theta = vector of hyperparameters: 
%       [\phi_R  \tau  \phi_{SS}  \Delta_{ss}  \phi_{S2S}  \Delta_{S2S}]
%       \phi_{R}: standard deviation of the model (Y)
%       \tau: standard deviation of B (event-effect)
%       \phi_{SS}: standard deviation of event-staion term 
%       \Delta_{SS} : decay parameter of event-station term
%       \phi_{S2S}: standard deviation of site term
%       \Delta_{S2S}  : decay parameter of site term
% NOTE:
%              can be matrix where each parameter is column vector of
%              length K, then output p is also column vector of length K 
%       lam_theta = prior lambda values:
%             [lam_phi_{R},lam_tau,lam_phi_{SS},lam_Delta_{SS},lam_phi_{S2S},lam_Delta_{S2S}]
%       dij = distance matrix (site to site)
%       W   = [X'*Z1 Z1 Z2]
%       mu_beta: vector of mean values of model coefficients (Beta)
%       var_beta: vector of standard deviation values of model coefficients
%       Note that Beta = [B1 B2 B3 B4] is defined weakly-informative 
%       zp: indx of observations
% OUTPUT p = scalar or vector (depends on input) of log-posterior
%            probability for given data and prior parameters
% -----------------------  NOTE & PROOF ** --------------------------------
% Extra explanation
% let CY = W*cov_eta*W' + cov_y
% (Get & Store) Cholesky decomposed matrix  : CCY  = chol(CY)
% for calculating the |CY| = prod(diag(CCY))^2 
% and the log(det(CY)) is 2*sum(log(diag(chol(CY))));
% computes YWM'/CY     is (CCY\(CCY'\YWM))'
% PROOF:
% A = G'*G, by definition of the Cholesky root
% log(det(A)) = log(det(G'*G)) = log(det(G')*det(G)) = 2*log(det(G))
% Since G is triangular, det(G) = prod(vecdiag(G))
% Therefore log(det(G))=sum(log(vecdiag(G)))
% Consequently, log(det(A)) = 2*sum(log(vecdiag(G)))


N  = numel(Y);         
NS = size(dij,1);     
NW = size(W,2);       
NP = length(mu_beta);
NT = NW-NS-NP;
K  = size(theta,1);

% hyperparameters
% Reparametrization using log(\theta) instead of \theta  
log_phi_R = theta(:,1);   
log_tau = theta(:,2); 
log_phi_SS = theta(:,3); 
log_Del_SS = theta(:,4); 
log_phi_S2S = theta(:,5); 
log_Del_S2S = theta(:,6); % Fixed for modeling 

% constant values for prior distribution ~ Eq.(A.12)
lam_Phi_R = lam_theta(1);
lam_tau = lam_theta(2);
lam_Phi_SS = lam_theta(3);
lam_Del_SS = lam_theta(4);
lam_Phi_S2S = lam_theta(5);
p = nan(K,1);

YWM   =  Y - W * [ mu_beta ; zeros(NT+NS,1)];  % used for (Y-W*MUe)

% Sigma_beta
cov_beta  = diag(var_beta);    

for a = 1:K
        % phi_{S2S} (Matern cov. of site terms) Eq.(3) & Eq.(A.8)
        cov_S2S = exp(2*log_phi_S2S(a))*exp(-dij/exp(log_Del_S2S(a)));
        
        % sub_phi_{SS} (Sub-Matern cov. of event-station terms) Eq.(4)
        cov_SS_sub = exp(2*log_phi_SS(a)).*(1+dij/exp(log_Del_SS(a))).*...
            exp(-dij/exp(log_Del_SS(a)));  
        
        % event_station (cov. of event-station terms) 
        cov_SS = kron(eye(NT,NT),cov_SS_sub); 
        cov_SS = cov_SS(zp,zp);
        % tau^2*I: covariance of log(\mu_es)
        cov_tau = exp(2*log_tau(a))*eye(NT);    
        % Phi^2_{R}*I : covariance of error data 
        cov_phi_R  = exp(2*log_phi_R(a))*eye(N);          
        % Sigma_eta: latent prior covariance Eq.(A.7), Eq.(A.16)
        cov_eta = [cov_beta      zeros(NP,NT)  zeros(NP,NS); 
                   zeros(NT,NP)  cov_tau       zeros(NT,NS);
                   zeros(NS,NP)  zeros(NS,NT)  cov_S2S  ]; 
        cov_y = cov_phi_R + cov_SS;  % Eq.(A.17)   
        %------------------------------------------------------------------
        CY = W*cov_eta*W' + cov_y;   % Eq.(A.19)
        CCY = chol(CY);    
        %------------------------------------------------------------------
        % log posterior of hyper parameters  Eq.(A.22)
        p(a) = -2*sum(log(diag(CCY)))/2 ...  
               - 0.5*((CCY\(CCY'\YWM))'*YWM) ...
               + log(lam_Phi_R)  - lam_Phi_R*exp(log_phi_R(a))   + log_phi_R(a) ...
               + log(lam_Phi_SS) - lam_Phi_SS*exp(log_phi_SS(a)) + log_phi_SS(a) ...
               + log(lam_Del_SS) - lam_Del_SS*exp(log_Del_SS(a)) + log_Del_SS(a) ...
               + log(lam_tau)    - lam_tau*exp(log_tau(a))       + log_tau(a)    ...
               + log(lam_Phi_S2S) - lam_Phi_S2S*exp(log_phi_S2S(a)) + log_phi_S2S(a); 
end
