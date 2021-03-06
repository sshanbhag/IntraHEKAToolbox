%-------------------------------------------------------------------------
% reading_data_example.m
%-------------------------------------------------------------------------
% IntraHEKAToolbox
%-------------------------------------------------------------------------
% code to demonstrate how to use the IntraHEKAToolbox to import
% intracellular data
% 
%-------------------------------------------------------------------------
% See Also: readPGF, readPUL, readDAT, readHEKA, HEKAview
%-------------------------------------------------------------------------

%-------------------------------------------------------------------------
% Sharad J. Shanbhag 
% sshanbhag@neomed.edu
%-------------------------------------------------------------------------
% Revisions: 
% 	2012 (SJS): Created 
% 	9 Jun 2016 (SJS): 
% 	 -	comments added
% 	 -	enabled git tracking
% 	 -	renamed to reading_data_example
%-------------------------------------------------------------------------

%-------------------------------------------------------------------------
% Define file names, paths to files
%-------------------------------------------------------------------------
% path to raw data directory
dpath = '/Users/sshanbhag/Work/Data/Mouse/IntracellularAmygdala/RawData/2012-03-12';
% data file name
dname = '2012-03-12-9.dat';

% pgf and pul files from HEKA software
%	pgf file has information about recording such as sample interval,
%	# sweeps, triggers, scaling, etc.
%
%	pul file has information about data series, such as resistance of
%	pipette, potential, etc.
pgfname = '2012-03-12-9.pgf';
pulname = '2012-03-12-9.pul';

% output file: dsc returned by readHEKA does not contain data - 
% instead, to conserve disk space, has a list of location of start 
% of each sweepwithin the data file (.dat) 
dscname = '/Users/sshanbhag/Desktop/2012-03-12-9.dsc';

%-------------------------------------------------------------------------
% Read in pgf information
%-------------------------------------------------------------------------
pgf = readPGF(fullfile(dpath, pgfname));

%-------------------------------------------------------------------------
% read pul file information
%-------------------------------------------------------------------------
[pul_tr, pul_rr] = readPUL(fullfile(dpath, pulname));

%-------------------------------------------------------------------------
% read data
%-------------------------------------------------------------------------
d = readDAT(fullfile(dpath, dname), pgf, pul_rr);

%-------------------------------------------------------------------------
% alternative method: readHEKA combines all steps
%-------------------------------------------------------------------------
[D, PGF, PUL, DSC] = readHEKA(fullfile(dpath, dname));

%-------------------------------------------------------------------------
% save data to file (mat file, but with .dsc extension)
%-------------------------------------------------------------------------
save(dscname, 'DSC', '-MAT');
