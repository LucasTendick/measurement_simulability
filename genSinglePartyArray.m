function SingleParty = genSinglePartyArray(oa,ma)
%GENSINGLEPARTYARRAY generates deterministic probability distributions
% This function has two required inputs
%   oa: the number of outcomes for each measurement
%   ma: the number of measurements
%
% SingleParty = genSinglePartyArray(oa,ma) generates a 3-D array,
% containing the ndet = oa^ma determinstic single-party probability
% distributions which can be generated when Alice has ma measurements, each
% with oa outcomes. P(a,x,lambda) = P(a|x,lambda) the probability to give
% as outcome a, when the measurement choice is x, given that the hidden
% variable has value lambda.
%
%   requires: nothing
%   authors: Paul Skrzypczyk, Daniel Cavalcanti
%   last updated: March 17, 2016

Ndet = oa^ma; %number of deterministic behaviours

SingleParty = zeros(oa,ma,Ndet); % initialise array

for lam = 0:Ndet-1
    lamdec = dec2base(lam,oa,ma)-'0'; % generates the string of outcomes a 
                                %(for each x), for the given variable lam
    for x = 1:ma
        for a = 0:oa-1
            SingleParty(a+1,x,lam+1) = (lamdec(x) == a);
            % probability = 1 if a = lamdec(x), 0 otherwise
        end
    end
end

% de2bi(lam,ma,oa);
end