function out = readLog(logfile, outfmt)
%-----------------------------------------------------------------------------
% out = readLog(logfile, outfmt)
%-----------------------------------------------------------------------------
% IntraHEKAToolbox
%-----------------------------------------------------------------------------
% Read log file
%
% Log file format (6 April, 2012):
% 
% 	comma-delimited file, 1 header line with columns:
% 	
% 		year
% 		month
% 		day
% 		unit_number
% 		animal
% 		cell_depth
% 		trial_start
% 		trial_end
% 		auditory_stim
% 		other_stim
% 		comments
% 
% 	
% 
%-----------------------------------------------------------------------------
% Input Arguments:
% 	logfile			data log file (usually .csv file)
%	logfmt			string for output format (default == 'struct'
% 							'struct'		struct array with singleton fields
% 							'array'		struct with cell array fields
% 
% Output Arguments:
%	out		[Ndatalines X 1] output struct array with fields 
%						determined from log file header line
%-----------------------------------------------------------------------------
% See also: 
%-----------------------------------------------------------------------------

%-----------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neomed.edu
%-----------------------------------------------------------------------------
% Created: 6 April, 2012 (SJS)
%
% Revisions:
%	11 Apr 2012 (SJS):
% 		- functionalized
% 		- added code to read in information and store in outstruct
%	24 Apr 2012 (SJS):
% 		- added output format selection
%-----------------------------------------------------------------------------
% TO DO:
%-----------------------------------------------------------------------------

% if no file provided, get one from user
if ~exist('logfile', 'var')
	[logfile, logpath] = uigetfile('*.csv', 'Select log file');
	if isequal(logfile, 0) || isequal(logpath, 0)
		outstruct = [];
		return
	end
	logfile = fullfile(logpath, logfile);
end

fprintf('Reading log file %s...\n', logfile);

% first, find how many lines there are in the log file
Nlines = countTextFileLines(logfile);
Ndatalines = Nlines - 1;
fprintf('\t%d total lines\n', Nlines);
fprintf('\t%d data lines\n\n', Ndatalines);

% open file
fp = fopen(logfile, 'rt');

% read fields and scan for fieldnames
tmpln = fgetl(fp);
tmp = textscan(tmpln, '%s', 'delimiter', ',');
fnames = tmp{1};
Nfields = length(fnames);

fprintf('\tFound %d fields:\n', Nfields);
for n = 1:Nfields
	fprintf('\t\t%s\n', fnames{n});
end
fprintf('\n\n');

% preallocate structure and string used to scan lines
txtscanstr = '';
for n = 1:Nfields
	tmpstruct.(fnames{n}) = [];
	txtscanstr = [txtscanstr '%s'];
end
outstruct = repmat(tmpstruct, Ndatalines, 1);

% loop through # of data lines
for n = 1:Ndatalines
	tmpln = fgetl(fp);
	tmp = textscan(tmpln, txtscanstr, 'delimiter', ',');
	% assign data to output struct array
	for f = 1:Nfields
		if ~isempty(tmp{f})
			outstruct(n).(fnames{f}) = tmp{f}{1};
		else
			outstruct(n).(fnames{f}) = [];
		end
	end
end

% close file
fclose(fp);

if ~exist('outfmt', 'var')
	outfmt = 'STRUCT';
end

switch upper(outfmt)
	
	case 'STRUCT'
		out = outstruct;
	case 'ARRAY'
		
		for f = 1:Nfields
			out.(fnames{f}) = cell(Ndatalines, 1);
			
			for n = 1:Ndatalines
				out.(fnames{f}){n} = outstruct(n).(fnames{f});
			end
		end

	otherwise
		warning('%s: unknown outfmt %s, using default (struct)', mfilename, outfmt);
		out = outstruct;
end




% ¡done!
