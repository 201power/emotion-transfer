function [label, model, llh] = emgms(X, init,m,n)
% EM algorithm for Gaussian mixture model
% Written by Mo Chen (mochen@ie.cuhk.edu.hk). March 2009.
% modifed by Li He 2011-3-7 to include spatial information
% for image segmentation used in color transfer
%% initialization
fprintf('EM Spatial for Gaussian mixture: running ... ');
R = initialization(X,init);
P=R;

tol = 1e-6;
maxiter = 500;
llh = -inf(1,maxiter);
converged = false;
t = 1;
while ~converged && t < maxiter
    t = t+1;
    model = maximization(X,R);
    [R1 llh(t)] = expectation(X,P,model,t);% LH
    P=bilateral(X,R,m,n,7);
    %s = sum(sum(R1-R));
    %converged = s < tol;
    converged = llh(t)-llh(t-1) < tol*abs(llh(t));
    R=R1;
end
[~,label(1,:)] = max(R,[],2);
llh = llh(2:t);
if converged
    fprintf('converged in %d steps.\n',t);
else
    fprintf('not converged in %d steps.\n',maxiter);
end

end

function R = initialization(X, init)
[d,n] = size(X);
if isstruct(init)  % initialize with model
    R  = expectation(X,init);
elseif length(init) == 1  % random initialization
    k = init;
    idx = randsample(n,k);
    m = X(:,idx);
    [~,label] = max(bsxfun(@minus,m'*X,sum(m.^2,1)'/2));
    while k ~= unique(label)
        idx = randsample(n,k);
        m = X(:,idx);
        [~,label] = max(bsxfun(@minus,m'*X,sum(m.^2,1)'/2));
    end
    R = full(sparse(1:n,label,1,n,k,n));
elseif size(init,1) == 1 && size(init,2) == n  % initialize with labels
    label = init;
    k = max(label);
    R = full(sparse(1:n,label,1,n,k,n));
elseif size(init,1) == d && size(init,2) > 1  %initialize with only centers
    k = size(init,2);
    m = init;
    [~,label] = max(bsxfun(@minus,m'*X,sum(m.^2,1)'/2));
    R = full(sparse(1:n,label,1,n,k,n));
else
    error('ERROR: init is not valid.');
end
end

function [R llh] = expectation(X, P, model,t)
mu = model.mu;
Sigma = model.Sigma;
w = model.weight;

n = size(X,2);
k = size(mu,2);
R = zeros(n,k);

for i = 1:k
    R(:,i) = loggausspdf(X,mu(:,i),Sigma(:,:,i));
end
R = bsxfun(@plus,R,log(w));
T = logsumexp(R,2);
llh = sum(T)/n; % loglikelihood
R = bsxfun(@minus,R,T);
R = exp(R);

if (t==2) P=R; end% first iteration
R = (R+P)/2;

end

function model = maximization(X, R)
[d,n] = size(X);
k = size(R,2);
sigma0 = eye(d)*(1e-6); % regularization factor for covariance

s = sum(R,1);
w = s/n;
mu = bsxfun(@rdivide, X*R, s);
Sigma = zeros(d,d,k);
for i = 1:k
    Xo = bsxfun(@minus,X,mu(:,i));
    Xo = bsxfun(@times,Xo,sqrt(R(:,i)'));
    Sigma(:,:,i) = (Xo*Xo'+sigma0)/s(i);
end

model.mu = mu;
model.Sigma = Sigma;
model.weight = w;
end