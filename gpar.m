function [upper,neff,lag1_corr,acpt_rate] = gpar(r)


% SOME FUNCTIONALITIES OF THE ORIGINAL FUNCTION HAVE BEEN REMOVED
% ONLY RELEVANT FUNCTIONALITIES FOR THE BHM2020 PROJECT WERE NOT REMOVED


% input :
% r = a n times m matrix mcmc iteration, that is m chains of length n

% output :
% upper = the Gelman--Rubin statistic
% neff  = the estimated effective sample size
% lag1_corr = the lag 1 autocorrelation
% acpt_rate = the acceptance ratio

% Written 1 July 1991 by Andrew Gelman, Dept. of Statistics, UC Berkeley.
% Function "monitor" added 28 Feb 1994, fixed 20 Jan 1995.
% (df+3)/(df+1) fix on 2 Oct 1995--see Brooks and Gelman (1995) paper.

% Please mail all comments/questions to gelman@stat.berkeley.edu

% Software can be freely used for non-commercial purposes and freely
%    distributed.  Submitted to statlib on 29 April 1992.  Current
%    version submitted to statlib on 29 Dec 1995.

% Derivations appear in the articles, "Inference from iterative simulation
%    using multiple sequences," by A. Gelman and D. B. Rubin, in "Statistical
%    Science" 7, 457-511 (1992), and "Some new ideas on inference from
%    iterative simulation using multiple sequences," by S. Brooks and
%    A. Gelman, technical report (1995).

%                     How to Use These Programs
%
%
% Preparation:  The results of m multiple sequences, each of length 2n,
%    of iterative simulation (e.g., the Metropolis algorithm or the
%    Gibbs sampler) used to simulate a probability distribution.
%    The starting points of the simulations should be drawn from a
%    distribution that is overdispersed relative to the "target"
%    distribution you are trying to simulate.
%
%    At each iteration, all scalar parameters or summaries of interest
%    should be stored.  The result will be an 2n by m matrix for each
%    scalar summary.
%
%    The program "gpar" (or "gpar.log" or "gpar.logit", see below) can be
%    used for each scalar summary, or the results from k scalar summaries
%    can be combined into an array of dimensions (2n, m, k) and put into
%    the program "monitor".
%
% To run:  Use gpar (or gpar.log or gpar.logit) to do output analysis
%    on the simulations for each scalar summary of interest.  For
%    each scalar summary r, gpar produces two results:
%
%       a.  Estimates of the distribution of r,
%       b.  The estimated potential scale reduction;
%           i.e., the factor by which the estimated t intervals for r
%           might decrease if the iterative simulation is continued
%           forever.  As n increases, the scale reduction approaches 1.

% The following S functions are included:

% gpar (Gibbs-parallel) takes a matrix r of simulations:  each of
%       m columns is an independent iteratively simulated sequence of length
%       2n.  The first n steps of each sequence are discarded.
%       Output is a list of three vectors:
%
%       post:  (2.5%, 50%, 97.5%) quantiles for the target distribution
%               based on the Student-t distribution
%       quantiles:  (2.5%, 25%, 50%, 75%, 97.5%) empirical quantiles of
%               the mn simulated values.
%       confshrink:  50% and 97.5% quantiles of a rough upper bound on
%               how much the confidence interval of "post" will shrink
%               if the iterative simulation is continued forever.
%
%           If both components of confshrink are not near 1, you should
%           probably run the iterative simulation further.


% gpar.log and gpar.logit are versions of gpar to be used for scalar
%       summaries that are all-positive (gpar.log) or are constrained to lie
%       between 0 and 1 (gpar.logit).

% monitor is a routine for monitoring the convergence of k scalar summaries
%       at once:  the input is an array of dimensions (2n, m, k).
%       Output is a k by 7 matrix, with 2.5%, 25%, 50%, 75%, 97.5% quantiles,
%       R-hat, and 97.5% upper bound for R-hat, for eack of the k scalar
%       summaries.
%       Optional additional input is trans:  a vector of length k, each
%       element of which is "", "log", or "logit", corresponding to no
%       transformation, log transformation, or logit transformation for
%       computing R-hat (the transformations have no effect on the quantiles).

%############################################################################

% MAIN PROGRAM

%function [confshrinkrange,postrange,quantiles]=gpar(r)
%function [upper,neff,lag1_corr,acpt_rate]=gpar(r)

%gpar _ function (r) {

if (max(size(size(r)))==2)
    alpha=0.05;                     % 95% intervals
    [n1,m1]=size(r);
    if (n1<m1),r=r';n=m1;m=n1;
        else n=n1;m=m1; 
    end        

% Compute the autocorrelation.
for mm=1:m,
    a0=corrcoef(r(2:n,mm),r(1:n-1,mm));
    lag1m(mm)=a0(1,2);
end
lag1_corr=mean(lag1m);

% Compute the acceptance ratio.
for mm=1:m,
    acpt(mm)=(length(find(diff(r(:,mm))))+1)/n;
end
acpt_rate=mean(acpt');

%        m _ ncol(r)
%        x _ r [(nrow(r)/2+1):nrow(r),]  # second half of simulated sequences
%        n _ nrow(x)
    x=r(floor(n1/2+1):n1,:);
    n=max(size(x));
        
% We compute the following statistics:
%
%  xdot:  vector of sequence means
%  s2:  vector of sequence sample variances (dividing by n-1)
%  W = mean(s2):  within MS
%  B = n*var(xdot):  between MS.
%  muhat = mean(xdot):  grand mean; unbiased under strong stationarity
%  varW = var(s2)/m:  estimated sampling var of W
%  varB = B^2 * 2/(m+1):  estimated sampling var of B
%  covWB = (n/m)*(cov(s2,xdot^2) - 2*muhat*cov(s^2,xdot)):
%                                               estimated sampling cov(W,B)
%  sig2hat = ((n-1)/n))*W + (1/n)*B:  estimate of sig2; unbiased under
%                                               strong stationarity
%  quantiles:  emipirical quantiles from last half of simulated sequences
    


    for jj=1:m
        xdot(jj)=mean(x(find(1-isnan(x(:,jj))),jj));  
        s2(jj)=var(x(find(1-isnan(x(:,jj))),jj));  
    end

    W=mean(s2);

    B=n*var(xdot);

    muhat=mean(xdot);

    varW=var(s2)/m;

    varB=B^(2)*2/(m-1);

    cov_s2_xdot_sqrd = cov(s2,xdot.^2);
    cov_s2_xdot = cov(s2,xdot);
    covWB=(n/m)*(cov_s2_xdot_sqrd(1,2)-2*muhat*cov_s2_xdot(1,2));

    sig2hat=((n-1)*W+B)/n;

    quantiles=prctile(x,100*[0.025,0.25,0.5,0.75,0.975]);      
   
    
    if (W > 1.e-8)   % non-degenerate case

        % Posterior interval post.range combines all uncertainties
        % in a t interval with center muhat, scale sqrt(postvar),
        % and postvar.df degrees of freedom.
        %
        %       postvar = sig2hat + B/(mn):  variance for the posterior interval
        %                               The B/(mn) term is there because of the
        %                               sampling variance of muhat.
        %       varpostvar:  estimated sampling variance of postvar


        postvar=sig2hat+B/(m*n);
        varpostvar= ...
        (((n-1)^2)*varW+(1+1/m)^2*varB+2*(n-1)*(1+1/m)*covWB)/n^2;

        postdf= 2*(postvar^2/varpostvar);
        postrange=muhat+sqrt(postvar)*tinv(1-alpha/2,postdf)*[-1,0,1];   

% Estimated potential scale reduction (that would be achieved by
% continuing simulations forever) has two components:  an estimate and
% an approx. 97.5% upper bound.
%
% confshrink = sqrt(postvar/W),
%     multiplied by sqrt(df/(df-2)) as an adjustment for the
%##      CHANGED TO sqrt((df+3)/(df+1))
%     width of the t-interval with df degrees of freedom.
%
% postvar/W = (n-1)/n + (1+1/m)(1/n)(B/W); we approximate the sampling dist.
% of (B/W) by an F distribution, with degrees of freedom estimated
% from the approximate chi-squared sampling dists for B and W.  (The
% F approximation assumes that the sampling dists of B and W are independent;
% if they are positively correlated, the approximation is conservative.)
        
         varlodf=2*(W^2/varW);


        confshrinkrange=sqrt([postvar/W, ...
                (n-1)/n+(1+1/m)*(1/n)*(B/W)*finv(0.975,m-1,varlodf)]* ...
                (postdf+3)/(postdf+1));
        upper=confshrinkrange(1,2);
        lower=confshrinkrange(1,1);

    end
elseif (max(size(size(r)))==3)
    alpha=0.05;                     % 95% intervals
    [ace,n1,m1]=size(r);
    if (n1<m1),r=r';,n=m1;,m=n1;
        else n=n1;,m=m1; 
    end        
    
    r=reshape(r(1,:,:),[n1 m]);
    % Compute the autocorrelation.
    for mm=1:m,
        a0=corrcoef(r(2:n,mm),r(1:n-1,mm));
        lag1m(mm)=a0(1,2);
    end
    lag1_corr=mean(lag1m);

   % Compute the acceptance ratio.
    for mm=1:m,
        acpt(mm)=(length(find(diff(r(:,mm))))+1)/n;
    end
    acpt_rate=mean(acpt');
    
    x=reshape(r((n1/2+1):n1,:),[length((n1/2+1):n1) m]);
    n=max(size(x));
    
    for jj=1:m
        xdot(jj)=mean(x(find(1-isnan(x(:,jj))),jj));  
        s2(jj)=var(x(find(1-isnan(x(:,jj))),jj));  
    end
    
    W=mean(s2)

    B=n*var(xdot);

    muhat=mean(xdot);

    varW=var(s2)/m;

    varB=B^(2)*2/(m-1);
    for j=1:m1,
        s20(j)=s2(1,j);
        xdot0(j)=xdot(1,j);
    end

    
    cov_s20_xdot0_sqrd = cov(s20,xdot0.^2);
    cov_s20_xdot0 = cov12(s20,xdot0);
    
    covWB=(n/m)*(cov_s20_xdot0_sqrd(1,2)-2*muhat*cov_s20_xdot0(1,2));

    sig2hat=((n-1)*W+B)/n;

   quantiles=prctile(x,100*[0.025,0.25,0.5,0.75,0.975]);      
%    quantiles=quantile(x,[0.025,0.25,0.5,0.75,0.975]);     

    if (W > 1.e-8)   % non-degenerate case

        postvar=sig2hat+B/(m*n)

        varpostvar= ...
            (((n-1)^2)*varW+(1+1/m)^2*varB+2*(n-1)*(1+1/m)*covWB)/n^2;

        postdf=2*(postvar^2/varpostvar);
                
        postrange=muhat+sqrt(postvar)*tinv(1-alpha/2,postdf)*[-1,0,1];   

        varlodf=2*(W^2/varW);
        
        confshrinkrange=sqrt([postvar/W, ...
            (n-1)/n+(1+1/m)*(1/n)*(B/W)*finv(0.975,m-1,varlodf)]* ...
        (postdf+3)/(postdf+1));
        upper=confshrinkrange(1,2);
        lower=confshrinkrange(1,1);
    end
end

% Compute the effective sample size.
neff=min(floor(n*m*postvar/B),n*m);


