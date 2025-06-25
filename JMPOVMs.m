function [are_JM_POVMs,Glam] = JMPOVMs(Max)
%JMPOVMS Determines whether a set of POVMs is jointly measurable or not
%  This function has one argument:
%   Max: a 4-D array, containing the POVM elements. The first two
%   dimensions contain the POVM elements, while the remaining two
%   dimensions are (a,x), such that Max(:,:,a,x) = M_a|x.
%
%  are_JM_POVMs = JMPOVMs(Max) is the indicator function for jointly
%  measurable POVMs. It returns 1 if the measurements are jointly
%  measurable , and 0 otherwise.
%
% [are_JM_POVMs,Glam] = JMPOVMs(Max) also returns in Glam the parent POVM
% which reproduces the POVMs Max. Glam is a dA x dA x ndet array, where
% ndet = oa^ma is the number outcomes of the parent POVM Glam,
% corresponding to all possible strings of outcomes for all ma
% measurements.
%
% are_JM_POVMs = JMPOVMs(Max) can also be used inside CVX as 
% a partially specified problem, to enforce the constraint that the CVX
% variable Max should be a set of jointly measurable POVMs
%
% EXAMPLE:
%   cvx_begin
%
%       variable Max(dA,dA,oa,ma)
%
%       subject to
%
%       JMPOVMs(Max) == 1
%
%   cvx_end
%
%  Inside CVX, Max is a set of measurements with ma inputs and oa outcomes
%  of dimension dA x dA. CVX enforces that Max should be jointly
%  measurable.
%
%   requires: CVX (http://cvxr.com/cvx/), QETLAB (http://www.qetlab.com)
%   authors: Paul Skrzypczyk, Daniel Cavalcanti
%   last updated: March 17, 2016


[dA,~,oa,ma] = size(Max);
% dA = dim., oa = # outcomes, ma = # inputs for Alice

Ndet = oa^ma; % number of deterministic strategies for Alice
SingleParty = genSinglePartyArray(oa,ma); % generate array containing 
                                            %single party strategies 

if isa(Max,'cvx') == 0 % if sigma isn't a CVX variable
    % check that the POVMs are valid
    if  validPOVMs(Max) == 0
        error('POVMs are not valid')
    end  
end

cvx_begin sdp quiet
    
    variable Glam(dA,dA,Ndet) hermitian semidefinite
    % Glam is the parent POVM

    subject to
    
    Max == squeeze(sum(repmat(Glam,[1,1,1,oa,ma])...
        .*permute(repmat(SingleParty,[1,1,1,dA,dA]),[4,5,3,1,2]),3));
    % M_a|x = sum_lambda D(a|x,lambda) G_lambda
    
    sum(Glam,3) == eye(dA);
    % sum_lambda G_lambda == 1

cvx_end

% CVX will return +inf if the problem is infeasible, and 0 if feasible
% this maps {+inf,0} to {0,1}
are_JM_POVMs = 1-min(cvx_optval,1);

% if there is no parent POVM, then return Glam as the empty array
if are_JM_POVMs == 0
    Glam = [];
end
    
end


