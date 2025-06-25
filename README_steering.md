### Code to accompany *[Quantum steering: a review with focus on semi-definite programming](https://doi.org/10.1088/1361-6633/80/2/024001)*
#### Daniel Cavalcanti and Paul Skrzypczyk

This repository provides a small collection of code which implements many of the semidefinite programs presented in the review article  
*[Quantum steering: a review with focus on semi-definite programming](https://doi.org/10.1088/1361-6633/80/2/024001)*      
Daniel Cavalcanti and Paul Skrzypczyk  
Rep. Prog. Phys. **80** 024001 (2017) .

All code is written in MATLAB and requires:
- [CVX](http://cvxr.com/) - a Matlab-based convex modeling framework
- [QETLAB](http://www.qetlab.com/) - A MATLAB Toolbox for Quantum Entanglement

It has been tested on Matlab R2014a, and CVX 2.1 

The code comprises the following:

- Membership, helper, and misc:
  - [NSAssemblage](https://github.com/paulskrzypczyk/steeringreview/blob/master/NSAssemblage.m): determine whether a bipartite assemblage is a valid non-signalling assemblage or not<sup>§</sup>.
  - [LHSAssemblage](https://github.com/paulskrzypczyk/steeringreview/blob/master/LHSAssemblage.m): determine whether a biparitte assemblage has an LHS model or not<sup>§</sup>.
  - [genAssemblage](https://github.com/paulskrzypczyk/steeringreview/blob/master/genAssemblage.m): generate an assemblage starting from a quantum state and a set of measurements.
  - [validPOVMs](https://github.com/paulskrzypczyk/steeringreview/blob/master/validPOVMs.m): determine whether a set of POVMs is valid or not<sup>§</sup>.
  - [JMPOVMs](https://github.com/paulskrzypczyk/steeringreview/blob/master/JMPOVMs.m): determine whether a set of measurements is jointly measurable or not<sup>§</sup>.
  - [genRandProjMeas](https://github.com/paulskrzypczyk/steeringreview/blob/master/genRandProjMeas.m): generate a random set of projective measurements.
  - [genSinglePartyArray](https://github.com/paulskrzypczyk/steeringreview/blob/master/genSinglePartyArray.m): generate the single-party determinstic probability distributions
  - [bestSteeringMeasurements](https://github.com/paulskrzypczyk/steeringreview/blob/master/bestSteeringMeasurements.m): find the optimal measurements given a state and a steering functional.
  - [bestSteeringState](https://github.com/paulskrzypczyk/steeringreview/blob/master/bestSteeringState.m): find the optimal state given a set of measurements and a steering functional.
  - [findRadiusPolytopeInBlochSphere](https://github.com/paulskrzypczyk/steeringreview/blob/master/findRadiusPolytopeInBlochSphere.m): determine the radius of the largest ball which can fit inside a polytope contained in the Bloch ball determined by a set of measurements<sup>¶</sup>.

- Steering quantifiers:
  - [steeringRobustness](https://github.com/paulskrzypczyk/steeringreview/blob/master/steeringRobustness.m): calculate the (standard/consistent) Steering Robustness of an assemblage.
  - [steeringWeight](https://github.com/paulskrzypczyk/steeringreview/blob/master/steeringWeight.m): calculate the (standard/consistent) Steering Weight of an assemblage.
  - [steeringLHSRobustness](https://github.com/paulskrzypczyk/steeringreview/blob/master/steeringLHSRobustness.m): calculate the (standard/consistent) Steering LHS-Robustness of an assemblage.
  - [steeringRobustnessState](https://github.com/paulskrzypczyk/steeringreview/blob/master/steeringRobustnessState.m): estimate the Steering Robustness of a state.
  - [steeringWeightState](https://github.com/paulskrzypczyk/steeringreview/blob/master/steeringWeightState.m): estimate the Steering Weight of a state.

- Local-Hidden-State models:
  - [targetStatePVMLHS](https://github.com/paulskrzypczyk/steeringreview/blob/master/targetStatePVMLHS.m): determine if a qubit-qudit state has an LHS model for all projective measurements on Alice.
  - [targetStatePOVMLHS](https://github.com/paulskrzypczyk/steeringreview/blob/master/targetStatePOVMLHS.m): determine if a qubit-qudit state has an LHS model for all POVMs on Alice.
  - [findPVMLHSStateGivenWitness](https://github.com/paulskrzypczyk/steeringreview/blob/master/findPVMLHSStateGivenWitness.m): find a qubit-qudit state that has an LHS model for all projective measurements on Alice and violates a given entanglement witness.
  - [findPOVMLHSStateGivenWitness](https://github.com/paulskrzypczyk/steeringreview/blob/master/findPOVMLHSStateGivenWitness.m): find a qubit-qudit state that has an LHS model for all POVMs on Alice and violates a given entanglement witness.  

- Multipartite steering:
  - [fLHSTripartite1Unt](https://github.com/paulskrzypczyk/steeringreview/blob/master/fLHSTripartite1Unt.m): determine if a tripartite assemblage with one untrusted device could have arisen from a fully-local state.
  - [fLHSTripartite2Unt](https://github.com/paulskrzypczyk/steeringreview/blob/master/fLHSTripartite2Unt.m): determine if a tripartite assemblage with two untrusted devices could have arisen from a fully-local state.
  - [bLHSTripartite1Unt](https://github.com/paulskrzypczyk/steeringreview/blob/master/bLHSTripartite1Unt.m): determine if a tripartite assemblage with one untrusted device could have arisen from a bi-separable state.
  - ~~bLHSTripartite2Unt~~: determine if a tripartite assemblage with two untrusted devices could have arisen from a bi-separable state. *still to come*

- Applications:
  - [localSteeringGuessProb](https://github.com/paulskrzypczyk/steeringreview/blob/master/localSteeringGuessProb.m): calculate the one-sided device-independent local guessing prob. of an assemblage.
  - [globalSteeringGuessProb](https://github.com/paulskrzypczyk/steeringreview/blob/master/globalSteeringGuessProb.m): calculate the one-sided device-independent global guessing prob. of an assemblage. 
  - [localSteeringGuessProbState](https://github.com/paulskrzypczyk/steeringreview/blob/master/localSteeringGuessProbState.m): estimate the one-sided device-independent local guessing prob. of a state.
  - [globalSteeringGuessProbState](https://github.com/paulskrzypczyk/steeringreview/blob/master/globalSteeringGuessProb.m): estimate the one-sided device-independent global guessing prob. of a state. 

<sup>§</sup>: These files can be used inside CVX as a means to enforce the corresponding constraint. 

<sup>¶</sup>: This file additionally needs [vert2lcon](http://www.mathworks.com/matlabcentral/fileexchange/30892).
