
% this function generates values with normal distribution 
% between interval [minVal maxVAl]
% Sahar rahpeyma
% 2014 
% University of Iceland

function x = randnLimit(mu, sigma, minVal, maxVal, varargin)

x = mu + sigma*randn(varargin{:});
% x = mu + sigma;
outsideRange = x<minVal | x>maxVal;
while nnz(outsideRange)>0
      a = find (outsideRange~=0); 
      for i= 1:length(a)
          x(a(i)) = mu(a(i)) + sigma.*randn(1,1);
%           x(a(i)) = mu(a(i)) + sigma;
          outsideRange = x<minVal | x>maxVal;
      end
end
