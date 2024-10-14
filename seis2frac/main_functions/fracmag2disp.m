function [dis, forward_mags, len]=fracmag2disp(mag,scaling_a,scaling_exp,mu,fracture_lengths)
%% FRACMAG2DISP
% Script to calculate displacement of a fracture from the seismic magnitude
% and fracture length, following Kanamori and Anderson (1975), where
% Displacement associated with fracturing comes from length-displacement
% relationships, where dmax = scaling_a * Length ^ scaling_exp
% (eg. dmax=0.0337L^1.02,
% https://www.sciencedirect.com/science/article/pii/S0191814119303657)
%
% Input:
%     mag: Seismic Magnitude of the fracturing
%     scaling_a: Constant, relating fracture length to Dmax
%     scaling_exp: Constant, relating fracture length to Dmax
%     mu: Shear Modulus
%     fracture_lengths: (Optional) Include known fracture lengths to
%                        calculate predicted magnitude using same law
%     
% Output:
%     dis: Displacement of fracture
%     forward_mags: Predicted magnitudes of fractures from given lengths
%     len: Length of fracture for each magnitude bin
%
% Jack McGrath July 21
    
mo=10.^((3/2)*(mag+6.07)); % Convert magnitude to moment
M=mo/(pi*mu*2*scaling_a); % Constant
len=nthroot(8*M,scaling_exp+2); % Rearrangement of Kanamori and Anderson equation to get length
dis=scaling_a.*(len.^scaling_exp); % Use scaling law to work out displacement

if nargin > 4
    forward_mags=(2/3)*log10((16/7)*(7/16)*pi*mu*(scaling_a*(fracture_lengths.^scaling_exp)/(fracture_lengths*0.5))*(0.5*fracture_lengths).^3)-6.07;
end
