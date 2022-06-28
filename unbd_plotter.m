% Ray Treinen, June 2022
%
% Compute the capillary surfaces known as symmetric unbounded liquid 
% bridges.  These surfaces have generating
% curves that contain vertial points at a point (sigma,T(sigma)) and are
% asypmtotic to the horizontal axis.
%
% This script needs Chebfun installed to run: chebfun.org
% The dependencies on chebfun are restricted to the generation of the
% spectral differentiation matrices and plotting.
%
% This script also needs the file "symmetric_capillary_unbounded" from this
% repository.  This will compute the data we need for computing the
% vertical points of the unbounded liquid bridges over an interval.
%
% If this code is used frequently, one can run the first line once and save
% the data in memory or a file, then the remaining lines can be run using
% that data.
%
% This code is slow because it is computing the function T(sigma) at 100
% points before computing the unbounded liquid bridge.  Consider the
% pervious paragraph if speed is of interest.

%% Physical parameters
% Choose the radius of the vertial point for the unbounded liquid bridge.
% The current range is from 0.085 to 2.0.  
% Larger values can be attained by altering the file 
% symmetric_capillary_unbounded
% but smaller values need a multi-scale algorithm.

Sigma = 1;

% Choose the terminal angle for the upper portion of the curve.  
% Note careully that phi0 is in [0,pi/2).

phi0 = 0;

%% Computationaly costly
% Uncomment for first running
% comment for faster performance while changing Sigma or phi0 on further
% runs of the code

[sig, Tvec, sig_min, sig_max ] = symmetric_capillary_unbounded();

%%
options = odeset('AbsTol',1e-11,'RelTol',1e-11);

Tsig = chebfun.interp1(sig,Tvec,'poly');
ic = [Sigma; Tsig(Sigma)];
[phi,vvv] = ode45(@FODE,[pi/2;phi0],ic,options);

figure(50)
plot(vvv(:,1),vvv(:,2),'k')
hold on

[phi,vvvv] = ode45(@FODE,[pi/2 .999*pi],ic,options);
plot(vvvv(:,1),vvvv(:,2),'k')
axis equal
