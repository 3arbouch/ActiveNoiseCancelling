function [f,g] = WeibAlphTAFCFitFun(x,beta,level,nCorrect,nError)
% [f,g] = WeibAlphTAFCFitFun(x,,beta,level,nCorrect,nError)
%
% 8/26/94		dhb		Wrote it.

% Unpack vector x
alpha = x(1);

% Force really big error if alpha or beta <= 0
if (alpha <= 0)
	f = 1e25;
	g = -1;
	return;
end

% Compute Weibull function on input levels
pCorrect = ComputeWeibTAFC(level,alpha,beta);

% Handle range problem, can't take log(0);
tol = 1e-4;
o_index = find(pCorrect == 1);
if (~isempty(o_index))
  pCorrect(o_index) = (1-tol)*ones(length(o_index),1);
end
  
% Compute error function, which is -log(likelihood).
tmp = nCorrect.*log(pCorrect) + nError.*log(1 - pCorrect);
f =  -sum(tmp);
g = -1;






