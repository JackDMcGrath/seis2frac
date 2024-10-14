function [MOC_ix,len_ix]=identify_fracture_sizes(eqs,bins,mag_window,MOC,len,vein_lengths)
%% IDENTIFY_FRACTURE_SIZES
%
% Script to identify the correct magnitude and lengths of fractures, and
% the bins that are above the magnitude of completeness

% Identify the bins of interest
mag_ix=find(((bins<=mag_window(2))+(bins>=mag_window(1)))==2); % Only look at earthquakes within the magnitude window
mag_ix=sort([mag_ix, mag_ix(1)-1, mag_ix(end)+1]);
MOC_ix=find(((bins<=max(eqs(:,4)))+(bins>=MOC))==2); % Get bins above MOC
MOC_ix=sort([MOC_ix, MOC_ix(end)+1]);
len_ix=find(((len<vein_lengths(2))+(len>vein_lengths(1)))==2); % Find the bins that are of the correct lengths