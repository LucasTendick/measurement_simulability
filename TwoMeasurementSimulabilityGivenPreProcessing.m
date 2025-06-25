function [GenRobTwoMeasurementSimulability, r] = TwoMeasurementSimulabilityGivenPreProcessing(Max,p)

[d,~,oa,ma] = size(Max);
% d = dim., oa = # outcomes, ma = # inputs for Alice

Ndet = oa^ma; % number of deterministic strategies for Alice
SingleParty = genSinglePartyArray(oa,ma); % create deterministic strategies for joint measurability test

% In the following, we use the characterization of 2-simulability according
% to Eq. (11) of the work 'Strict hierarchy between n-wise measurement simulability,
% compatibility structures, and multi-copy compatibility' (see: arXiv:XXXX.YYYYY)


cvx_begin 

variable Fax_1(d,d,oa,ma) hermitian semidefinite % The first jointly measuarable assemblage
variable Fax_2(d,d,oa,ma) hermitian semidefinite % The second jointly measuarable assemblage
variable Nax(d,d,oa,ma) hermitian semidefinite   % The noise model assemblage
variable r(1)

maximize r % The noise parameter

subject to

r >= 0; 

JMPOVMs(Fax_1) == 1; % Make sure it is JM
JMPOVMs(Fax_2) == 1; % Make sure it is JM
squeeze(sum(Nax,3)) == repmat((1-r)*eye(d),[1,1,ma]); % Normalization

for x = 1:ma
    for a = 1:oa
         Nax(:,:,a,x) == (1-r)*trace(Max(:,:,a,x))*eye(d)/(d);  
         % optional if we want the depolarized robustness, otherwise we
         % compute the generalized robustness
         r*Max(:,:,a,x)+Nax(:,:,a,x) == p(1,x)*Fax_1(:,:,a,x)+p(2,x)*Fax_2(:,:,a,x); % see Eq.(11)
    end
end



cvx_end


GenRobTwoMeasurementSimulability = (1/r)-1; % Conversion to usual representation of Generalized robustness. 


end