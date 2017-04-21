function [data, pgf, pul, dsc] = readHEKA(datfile, varargin)
%-----------------------------------------------------------------------------
% [data, pgf, pul] = readHEKA(datfile, varargin)
%-----------------------------------------------------------------------------
% 
% Reads data from HEKA .dat file given by datfile
%
% Optional arguments:
% 	'Downsample', <downsample_factor>
% 		
% 		downsamples data sweeps by factor of downsample_factor
% 		e.g., readHEKA('test.dat', 'Downsample', 2)
%
% Essentially, a wrapper around three functions:
% 		readPGF()
% 		readPUL()
% 		readDAT()
%		buildDSC()
% 
%-----------------------------------------------------------------------------
% Input Arguments:
% 	datfile		path/filename
% 	
% 	'Downsample', <factor>		downsamples data traces
%
% Output Arguments:
% 	data		trace data structure
%	pgf		pgf structure (from readPGF)
% 	pul		pul data (from readPUL)
% 		pul.tr
% 		pul.rr
%-----------------------------------------------------------------------------
% See also: readPGF, readPUL, readDAT
%-----------------------------------------------------------------------------

%-----------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neomed.edu
%-----------------------------------------------------------------------------
% Created: 22 March, 2012 (SJS)
% 				based on Visual Basic code clsOData.vb by Olga Galazyuk
%
% Revisions:
%	-	added downsample input parameter
%-----------------------------------------------------------------------------
% TO DO:
%-----------------------------------------------------------------------------

%-----------------------------------------------------------------------------
%-----------------------------------------------------------------------------
% Read information from file
%-----------------------------------------------------------------------------
%-----------------------------------------------------------------------------

if ~exist('datfile', 'var')
	% open gui dialog
	[filename, pathname] = uigetfile('*.dat', 'Pick a HEKA data file');
	if isequal(filename, 0) || isequal(pathname, 0)
		disp('Cancelled')
		data = [];
		pgf = [];
		pul = [];
		dsc = [];
		return
	else
		datfile = fullfile(pathname, filename);
		clear pathname filename
	end	
end

[dpath, dfile_base, ~] = fileparts(datfile);

datname = [dfile_base '.dat'];
pgfname = [dfile_base '.pgf'];
pulname = [dfile_base '.pul'];
dscname = [dfile_base '.dsc']; %#ok<NASGU>

% check files
if ~exist(fullfile(dpath, pgfname), 'file')
	error('%s: pgf file %s not found', mfilename, pgfname);
end
if ~exist(fullfile(dpath, pulname), 'file')
	error('%s: pul file %s not found', mfilename, pulname);
end
if ~exist(fullfile(dpath, datname), 'file')
	error('%s: dat file %s not found', mfilename, datname);
end

% read in pgf file
pgf = readPGF(fullfile(dpath, pgfname));

% read in pul file
pul = struct('tr', [], 'rr', []);
[pul.tr, pul.rr] = readPUL(fullfile(dpath, pulname));

% read in dat file
data = readDAT(fullfile(dpath, datname), pgf, pul.rr);

% build dsc struct
dsc = buildDSC(fullfile(dpath, datname), pgf, pul.rr);

% decimate (downsample) raw trace data if so instructed and 
% update sampling rate to reflect this
if nargin > 1
	vcheck = strcmpi('Downsample', varargin);
	if (any(vcheck) && length(varargin) > 1)
		vindx = find(vcheck) + 1;
		decifactor = varargin{vindx};
		
		for d = 1:length(data)
			data(d).trace1 = decimate(data(d).trace1, decifactor);
			data(d).trace2 = decimate(data(d).trace2, decifactor);
			data(d).rate = data(d).rate * decifactor;
		end
	end
end

