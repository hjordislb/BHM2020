function H = my_hessian(func,loc,h,varargin)
% FUNCTION : 
% This function calculate Hessian matrix used for BHM modeling by central
% finite differences any given function handle "func", which accepts parameters,
% any number of parameters (based on input "loc" and "h")
% *********************************************************************** %
% INPUT :
% func = function handle for function to evaluate (e.g. BHM_lnpost)
% loc  = vector for parameter location where Hessian is evaluated
% h    = step width vector (one element per parameter)
% varargin = additional data, if required by function
% *********************************************************************** %
% OUTPUT :
% H = Hessian Matrix
% *********************************************************************** %
% NOTE: input function must be of form: func(loc,varargin{:})
% *********************************************************************** %
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
%  DATE    : Nov. 2016               
% *********************************************************************** %

M = length(loc);
H = nan(M);
e = [1 1];
FX = func(loc,varargin{:}); % n=1
for a = 1:M    % 2nd order diff. w.r.t. same parameter
    H(a,a) = (func(loc + sparse(1,a,1,1,M).*h,varargin{:}) -2*FX ...
            + func(loc - sparse(1,a,1,1,M).*h,varargin{:}))/(h(a)^2);
end 
if M<2
    return
end
for a = 1:M-1     % diff. w.r.t. parameter a
    for b = a+1:M % diff. w.r.t. parameter b
        H(a,b) = (func(loc + sparse(e,[a b],e,1,M).*h,varargin{:})...
            -func(loc + sparse(e,[a b],[1 -1],1,M).*h,varargin{:})...
            -func(loc + sparse(e,[a b],[-1 1],1,M).*h,varargin{:})...
            +func(loc + sparse(e,[a b],-[1 1],1,M).*h,varargin{:}))...
            /(4*h(a)*h(b));
        H(b,a) = H(a,b);
    end
end 

