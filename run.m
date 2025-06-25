% This function is based on the work 'Strict hierarchy between n-wise measurement simulability,
% compatibility structures, and multi-copy compatibility' (see: arXiv:XXXX.YYYYY)


% This manuscript assists the proof of Result 1 and provides the necessary
% data of the 2-simulability with deterministic pre-processing and also
% produces the data corresponding to the (heuristically optimal)
% probabilistic 2-simulability strategy.


% Requires: CVX (http://cvxr.com/cvx/), A solver like SDPT3, Sedumi or
% Mosek, etc., Steeringreview-master (https://arxiv.org/abs/1604.00501) (https://github.com/paulskrzypczyk/steeringreview) 
% for the function genSinglePartyArray and JMPOVMs.
% Furthermore, it requires QETLAB: (https://github.com/nathanieljohnston/QETLAB, https://qetlab.com/)

% Authors: Lucas Tendick, Costantino Budroni, Marco Túlio Quintino
% Code by: Lucas Tendick
% Last Update: 25.06.2025

% The main function that is used here is the Function:
% TwoMeasurementSimulabilityGivenPreProcessing.


% Let us define the required measurements:

H = (1/sqrt(2))*[1 1; 1 -1];

Max(:,:,1,1) = (1/2)*(eye(2)+Pauli(1));
Max(:,:,2,1) = (1/2)*(eye(2)-Pauli(1));

Max(:,:,1,2) = (1/2)*(eye(2)+Pauli(3));
Max(:,:,2,2) = (1/2)*(eye(2)-Pauli(3));

Max(:,:,1,3) = (1/2)*(eye(2)+H);
Max(:,:,2,3) = (1/2)*(eye(2)-H);

% We first compute the visibility according to the probabilistic
% pre-processing

p = [1 0 1/2; 0 1 1/2];

[~,eta_prob] = TwoMeasurementSimulabilityGivenPreProcessing(Max,p);

% Now, let us handle the deterministic pre-processings. 
% We do this in the most naive way by simply listing all strategies and
% going through them.

q(:,:,1) = [1 1 1; 0 0 0];
q(:,:,2) = [1 1 0; 0 0 1];
q(:,:,3) = [1 0 1; 0 1 0];
q(:,:,4) = [1 0 0; 0 1 1];
q(:,:,5) = [0 1 1; 1 0 0];
q(:,:,6) = [0 1 0; 1 0 1];
q(:,:,7) = [0 0 1; 1 1 0];
q(:,:,8) = [0 0 0; 1 1 1];

for cnt = 1: size(q,3)

[~,r(cnt)] = TwoMeasurementSimulabilityGivenPreProcessing(Max,q(:,:,cnt));

end

eta_determ = max(r);
