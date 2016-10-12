%-----------------------------------------------------------------------------
% experiment
%-----------------------------------------------------------------------------
% IntraHEKA Toolbox
% Class Definition
%-----------------------------------------------------------------------------
% 
%-----------------------------------------------------------------------------
% See also: sweep (class), readPGF, readPUL, readHEKA, readDAT
%-----------------------------------------------------------------------------

%-----------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neomed.edu
%-----------------------------------------------------------------------------
% Created: 12 April, 2012 (SJS)
%
% Revisions:
%	19 Apr 2012 (SJS):
% 	 -	added PlotOpts property
%			has Color and Scale fields
% 	 -	moved plotting methods into separate files
%	 -	added something...
%	4 May 2012 (SJS):
% 	 -	added Statistics property to store statistic objects
% 
%-----------------------------------------------------------------------------
% TO DO:
%-----------------------------------------------------------------------------

classdef experiment < handle
	properties
		Sweeps
		Nsweeps
		PGF
		PUL
		BaseName
		BasePath
		DATfilename
		PGFfilename
		PULfilename
		EXPfilename
		Info
		PlotOpts
		isInitialized = 0;
		Statistics
		Nstatistics = 0;
	end
	
	methods
		
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		% Constructor
		%------------------------------------------------------------------------
		function obj = experiment(base_path, base_name, cmdstr)
			
			% default PlotOpts
			obj.PlotOpts.Color = 'b';
			obj.PlotOpts.Scale = 10;
			
			% if no arguments given, return an "empty" instance for most 
			% properties/values
			if nargin == 0
 				obj.Sweeps = [];
				obj.Nsweeps = [];
				obj.PGF = [];
				obj.PUL = [];
				obj.BaseName = '';
				obj.BasePath = '';
				obj.DATfilename = '';
				obj.PGFfilename = '';
				obj.PULfilename = '';
				obj.EXPfilename = '';
				obj.Info = [];
 				obj.Statistics = statistic;
				return
			end
			
			% little trick to get path and base name in canonical form
			[base_path, base_name, ext] = fileparts(fullfile(base_path, base_name)); %#ok<ASGLU>
			% otherwise, construct file names...
			obj.Sweeps = sweep;
			obj.Nsweeps = [];
			obj.PGF = [];
			obj.PUL = [];
			obj.BaseName = base_name;
			obj.BasePath = base_path;
			obj.DATfilename = [base_name '.dat'];
			obj.PGFfilename = [base_name '.pgf'];
			obj.PULfilename = [base_name '.pul'];
			obj.EXPfilename = [base_name '.exp'];
			obj.Info = [];
 			obj.Statistics = statistic;
			
			% ...and read in data if so instructed
			if nargin == 3
				switch upper(cmdstr)
					case 'INITIALIZE'
						obj.Initialize;
					otherwise
						warning('%s: unknown command %s', mfilename, upper(cmdstr));
				end
			end
		end
		%------------------------------------------------------------------------

		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		function Initialize(obj)
		%------------------------------------------------------------------------
		% initialize object (read in data)
		% must be called before doing most of the useful things in the 
		% experiment object.  
		%------------------------------------------------------------------------
			
			% check if BaseName is set/defined
			if isempty(obj.BaseName)
				% if not, throw warning and return
				warning('%s: no base_name defined!', mfilename)
				return
			end
			
			% make local versions of full file names
			datfile = fullfile(obj.BasePath, obj.DATfilename);
			pgffile = fullfile(obj.BasePath, obj.PGFfilename);
			pulfile = fullfile(obj.BasePath, obj.PULfilename);
			% check to make sure they can be found
			if ~exist(datfile, 'file')
				warning('%s: data file %s not found', mfilename, datfile);
				return
			elseif ~exist(pgffile, 'file')
				warning('%s: pgf file %s not found', mfilename, pgffile);
				return
			elseif  ~exist(pulfile, 'file')
				warning('%s: pul file %s not found', mfilename, pulfile);
				return
			end
			
			%----------------------------------------------------
			% read in various datums
			%----------------------------------------------------
			% read info from pgf file
			obj.PGF = readPGF(pgffile);
			% read info from pul file
			obj.PUL = struct('tr', [], 'rr', []);
			[obj.PUL.tr, obj.PUL.rr] = readPUL(pulfile);
				
			%----------------------------------------------------
			% build sweeps
			%----------------------------------------------------
			obj.BuildSweeps;
			
			% set the isInitialized property to 1 (true)
			obj.isInitialized = 1;
		end
		%------------------------------------------------------------------------
		
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		function BuildSweeps(obj)
		%------------------------------------------------------------------------
		% builds the Sweeps object array property.  
		% This is usually called by the experiment.Initialize() method, but
		% can be called independently if you know what you're doing...
		%------------------------------------------------------------------------
		
			% use buildDSC function to get sweep information
			dsc = buildDSC(fullfile(obj.BasePath, obj.DATfilename), obj.PGF, obj.PUL.rr);
			
			obj.Nsweeps = dsc.nsweeps;
			obj.Info = struct('Scale', [], 'Gain', []);
			obj.Info.Scale = [dsc.dfactor1, dsc.dfactor2];
			
			%%%%% Set Gain for channels
			%%%%% note that this is arbitrary and is done to get the
			%%%%% units to work out in Volts!!!!!!!!!!!!
			obj.Info.Gain = [1 10];
			
			% allocate Sweeps object array
			try
				obj.Sweeps = repmat(sweep, obj.Nsweeps, 1);
			catch
				for n = 1:obj.Nsweeps
					obj.Sweeps(n) = sweep;
				end
			end
			
			% loop through # of sweeps, assign values to sweeps
			for n = 1:obj.Nsweeps
				obj.Sweeps(n).Nsamples = dsc.size(n);
				obj.Sweeps(n).Sweepindex = n;
				obj.Sweeps(n).Rate = dsc.rate(n);
				obj.Sweeps(n).fp1 = dsc.trace1_fp(n);
				obj.Sweeps(n).fp2 = dsc.trace2_fp(n);
			end
		end
		%------------------------------------------------------------------------
		
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		function [trace1, trace2] = GetRawTrace(obj, sweeplist)
		%------------------------------------------------------------------------
		% get raw (unscaled) trace(s)
		% returns [trace1, trace2], where trace1 is usually the stimulus, 
		% and trace2 is the electrode data
		%------------------------------------------------------------------------
			if ~obj.isInitialized
				warning('%s: object not initialized', mfilename);
				return
			end
			if any(~between(sweeplist, 1, obj.Nsweeps))
				error('%s: sweeplist out of bounds (Nsweeps = %d)', mfilename, obj.Nsweeps);
			end
			datfile = fullfile(obj.BasePath, obj.DATfilename);
			ns = length(sweeplist);
			if ns == 1
				[trace1, trace2] = obj.Sweeps(sweeplist).ReadData(datfile);
				return
			end
			trace1 = cell(ns, 1);
			trace2 = cell(ns, 1);
			for n = 1:length(sweeplist)
				[trace1{n}, trace2{n}] = obj.Sweeps(sweeplist(n)).ReadData(datfile);
				trace1{n} = trace1{n};
				trace2{n} = trace2{n};
			end
		end
		%------------------------------------------------------------------------
		
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		function [trace1, trace2] = GetTrace(obj, sweeplist)
		%------------------------------------------------------------------------
		% get scaled trace(s)
		% returns [trace1, trace2], where trace1 is usually the stimulus, 
		% and trace2 is the electrode data
		%------------------------------------------------------------------------
			if ~obj.isInitialized
				warning('%s: object not initialized', mfilename);
				return
			end
			
			if any(~between(sweeplist, 1, obj.Nsweeps))
				error('%s: sweeplist out of bounds (Nsweeps = %d)', mfilename, obj.Nsweeps);
			end
			
			datfile = fullfile(obj.BasePath, obj.DATfilename);
			
			ns = length(sweeplist);

			if ns == 1
				[trace1, trace2] = obj.Sweeps(sweeplist).ReadData(datfile);
				trace1 = obj.Info.Scale(1) .* trace1 ./ obj.Info.Gain(1);
				trace2 = obj.Info.Scale(2) .* trace2 ./ obj.Info.Gain(2);
				return
			end
			trace1 = cell(ns, 1);
			trace2 = cell(ns, 1);
			for n = 1:length(sweeplist)
				[trace1{n}, trace2{n}] = obj.Sweeps(sweeplist(n)).ReadData(datfile);
				trace1{n} = obj.Info.Scale(1) .* trace1{n} ./ obj.Info.Gain(1);
				trace2{n} = obj.Info.Scale(2) .* trace2{n} ./ obj.Info.Gain(2);
			end
		end
		%------------------------------------------------------------------------
		
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		function [rate, trace1, trace2] = GetResampledTrace(obj, sweeplist, decifactor)
		%------------------------------------------------------------------------
		% get scaled and resampled trace(s)
		% returns [rate, trace1, trace2], where trace1 is usually the stimulus, 
		% and trace2 is the electrode data. rate is new sampling rate
		%------------------------------------------------------------------------
			% check decimation factor
			if decifactor < 1
				error('%s: decifactor must be integer >= 1 (decifactor = %d)', mfilename, decifactor);
			end
			% get the trace(s) in raw format
			[trace1, trace2] = obj.GetTrace(sweeplist);
			% if decifactor == 1, nothing to do!
			if decifactor == 1
				rate = obj.Sweeps(sweeplist(1)).Rate;
				return
			end
			% # of sweeps
			ns = length(sweeplist);
			% determine original sample rate - this is a bit of kludge, under
			% the possibility that sample rate might vary across sweeps, 
			% which is unlikely...  nevertheless...
			tmprate = zeros(size(sweeplist));
			for n = 1:ns
				tmprate(n) = obj.Sweeps(sweeplist(n)).Rate;
			end
			% use the lowest value of the found rates (actually, sample intervals,
			% in seconds) and multiply by decifactor to get the new rate
			rate = min(tmprate) * decifactor;
			% do the actual decimation
			if ns == 1
				trace1 = decimate(trace1, decifactor);
				trace2 = decimate(trace2, decifactor);
			else
				for n = 1:ns
					trace1{n} = decimate(trace1{n}, decifactor);
					trace2{n} = decimate(trace2{n}, decifactor);
				end
			end
		end
		%------------------------------------------------------------------------
		
		
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		function SetInfoFromLog(obj, logData)
		%------------------------------------------------------------------------
		% set stimulus information
		% apply log info to sweeps
		%------------------------------------------------------------------------
			% make sure object is initialized
			if ~obj.isInitialized
				warning('%s: object not initialized', mfilename);
				return
			end
			% first, need to go through and combine the year, month, day unit fields
			nLog = length(logData);
			unitNames = cell(nLog, 1);
			for l = 1:nLog
				L = logData(l);
				unitNames{l} = sprintf('%s-%s-%s-%s', L.year, L.month, L.day, L.unit_number);
			end
			% then use this to find matches
			logMatches = find(strcmpi(obj.BaseName, unitNames));
			nMatches = length(logMatches);
			if nMatches == 0
				% if no matches found, warn user
				warning('nMatches == 0!');
				return
			else
				% feedback on cmd line
				fprintf('Matches in log file\n');
				fprintf('\t%d\n', logMatches);
				% set animal and depth information
				obj.Info.Animal = logData(logMatches(1)).animal;
				obj.Info.Depth = str2num(logData(logMatches(1)).cell_depth); %#ok<*ST2NM>
				obj.Info.Nconditions = nMatches;
				% then, build Stimulus information struct array
				Stimulus = repmat(	struct(	'TrialStart', [], ...
														'TrialEnd', [], ...
														'Ntrials', [], ...
														'AuditoryStimulus', '', ...
														'OtherStimulus', '', ...
														'Comments', ''), ... 
											nMatches, 1);
				for n = 1:nMatches
					L = logData(logMatches(n));
					Stimulus(n).TrialStart = str2num(L.trial_start);
					Stimulus(n).TrialEnd = str2num(L.trial_end);
					Stimulus(n).Ntrials = 1 + Stimulus(n).TrialEnd - Stimulus(n).TrialStart;
					Stimulus(n).AuditoryStimulus = L.auditory_stim;
					Stimulus(n).OtherStimulus = L.other_stim;
					Stimulus(n).Comments = L.comments;
				end
				% store Stimulus in object
				obj.Info.Stimulus = Stimulus;
			end
			
		end
		%------------------------------------------------------------------------
		
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		function [t2avg, t2std, t1] = MeanTraceForCondition(obj, Condition, varargin)
		%------------------------------------------------------------------------
		% returns average trace for specific Condition
		%------------------------------------------------------------------------
			% make sure object is initialized
			if ~obj.CheckInitAndCondition(Condition)
				error('%s: Condition %d not found or class is not initialized', ...
								mfilename, Condition);
			end
			% check input args
			DESPIKE = 0;
			optargin = size(varargin,2);
% 			stdargin = nargin - optargin;
			if optargin
				% parse varargin args
				n = 1;
				while n < optargin
					if ischar(varargin{n})
						if strcmpi(varargin{n}, 'NoSpikes')
							DESPIKE = 1;
							n = n + 1;
						end
					end
				end
			end
			% get the sweeps for the indicated Condition
			sinfo = obj.GetStimulusForCondition(Condition);		
			sweeplist = obj.GetSweepListForCondition(Condition);
			[trace1, trace2] = obj.GetTrace(sweeplist);
			% store first trace1 data, then clear trace1 to save space
			t1 = trace1{1};
			clear trace1
			% despike traces if needed
			if DESPIKE
				for s = 1:length(sweeplist)
					trace2{s} = deSpike(trace2{s}, obj.Sweeps(sweeplist(s)).Rate);
				end
			end
			% compute mean of trace2 cell array
			t2avg = mean(cell2mat(trace2'), 2);
			t2std = std(cell2mat(trace2'), 0, 2);
		end
		%------------------------------------------------------------------------
		
		
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------		
		function tlen = GetTraceLengthForCondition(obj, Condition) 
		%------------------------------------------------------------------------
		% Returns trace lengths for a given condition
		%------------------------------------------------------------------------
			% make sure object is initialized
			if ~obj.CheckInitAndCondition(Condition)
				tlen = [];
				return
			end
			% get the sweeps for the indicated Condition
			sinfo = obj.GetStimulusForCondition(Condition);		
			sweeplist = obj.GetSweepListForCondition(Condition);
			tlen = zeros(length(sweeplist), 1);
			% get # of samples for each sweep
			for t = 1:length(tlen)
				tlen(t) = obj.Sweeps(sweeplist(t)).Nsamples;
			end
		end
		%------------------------------------------------------------------------
		
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		function sweeplist = GetSweepListForCondition(obj, Condition) 
		%------------------------------------------------------------------------
		% returns sweep list (indices) for given Condition
		%------------------------------------------------------------------------
			% make sure object is initialized
			if ~obj.CheckInitAndCondition(Condition)
				sweeplist = [];
				return
			end
			% get the sweeps for the indicated Condition
			sinfo = obj.GetStimulusForCondition(Condition);		
			sweeplist = sinfo.TrialStart : sinfo.TrialEnd;
		end
		%------------------------------------------------------------------------
		
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		function s = GetStimulusForCondition(obj, Condition)
		%------------------------------------------------------------------------
		% Method to get Stimulus information for specific condition
		%------------------------------------------------------------------------
			% make sure object is initialized
			if ~obj.CheckInitAndCondition(Condition)
				s = [];
				return
			else
				s = obj.Info.Stimulus(Condition);
			end
		end
		%------------------------------------------------------------------------

		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		function out = CheckInitAndCondition(obj, Condition)
		%------------------------------------------------------------------------
		% internal method to check initialization status and validity of 
		% Condition value
		%------------------------------------------------------------------------
			out = 0;
			% make sure object is initialized
			if ~obj.isInitialized
				warning('%s: object not initialized', mfilename);
				return
			elseif Condition == obj.Info.Nconditions
				out = 1;
				return
			% check that Condition is inbounds
			elseif ~between(Condition, 1, obj.Info.Nconditions)
				warning('%s: Condition %d not within bounds (Nconditions: %d)', ...
								mfilename, Condition, obj.Info.Nconditions);
				return
			else
				% all okay!
				out = 1;
				return
			end
		end
		%------------------------------------------------------------------------
		
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------		
		function stim = GetStimParamForCondition(obj, Condition, varargin) 
		%------------------------------------------------------------------------
		% Returns stimulus information for a given condition
		% Note that these are parameters that are computed directly from the 
		% stimulus trace and may be inaccurate!!!
		%------------------------------------------------------------------------
		%	Input Args:
		% 		Condition		Experiment Condition ID #
		% 
		% 		decifactor			decimation factor (default = 10)
		% 		threshold			threshold for start/end of stimulus
		% 								detections (default = 0.1);
		% 	Output Args:
		% 
		% 		stim	structure with fields:
		%				start		stimulus start time (milliseconds)
		%				end		stimulus end time (milliseconds)
		%------------------------------------------------------------------------
			% decimation factor
			DFACT = 10;
			% threshold for finding start/end
			THRESHOLD = 0.1;
			
			% make sure object is initialized
			if ~obj.CheckInitAndCondition(Condition)
				stim = [];
				return
			end
			
			% check inputs
			switch length(varargin)
				case 1
					if ~isempty(varargin{1})
						DFACT = varargin{1};
					end
				case 2
					if ~isempty(varargin{1})
						DFACT = varargin{1};
					end
					if ~isempty(varargin{2})
						THRESHOLD = varargin{2};
					end
			end

			% get the first stim trace for the indicated Condition
			sinfo = obj.GetStimulusForCondition(Condition);
			sweeplist = obj.GetSweepListForCondition(Condition);
			[rate, trace1, trace2] = obj.GetResampledTrace(sweeplist(1), DFACT);
			
			% then, compute envelope using Hilbert transform and find start and end
			env = abs(hilbert(1000*trace1));
			
			stim.start = 1000 * rate * find(env > THRESHOLD, 1, 'first');
			stim.end = 1000 * rate * find(env > THRESHOLD, 1, 'last');

		end
		%------------------------------------------------------------------------
		
		
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------		
		function [statIndex, statObj] = AddStatistic(obj, Statname, Statvals)
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		%	Input Args:
		% 	Output Args:
		%------------------------------------------------------------------------
			if obj.Nstatistics == 0
				obj.Nstatistics = 1;
			else
				obj.Nstatistics = obj.Nstatistics + 1;
			end
			obj.Statistics(obj.Nstatistics) = statistic(Statname, Statvals);
			statIndex = obj.Nstatistics;
			statObj = obj.Statistics(obj.Nstatistics);
			return
		end
		%------------------------------------------------------------------------
		
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------		
		function [statIndex, statObj] = GetStatisticByName(obj, Statname, CompareMethod)
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		%	Input Args:
		% 	Output Args:
		%------------------------------------------------------------------------
			if obj.Nstatistics == 0
				statIndex = 0;
				statObj = [];
				return
			end			

			if ~exist('CompareMethod', 'var')
				CompareMethod = 'STRCOMPI';
			end
			
			statnamelist = obj.GetStatisticNameList;
			
			switch upper(CompareMethod)
				case 'STRCMP'
					namematches = strcmp(Statname, statnamelist);
				case 'STRCMPI'
					namematches = strcmpi(Statname, statnamelist);
				case 'STRNCMP'
					namematches = strncmp(Statname, statnamelist, length(Statname));
				case 'STRNCMPI'
					namematches = strncmpi(Statname, statnamelist, length(Statname));
				otherwise
					error('%s: unknown stat string compare method %s', mfilename, CompareMethod);
			end
			
			if any(namematches)
				statIndex = find(namematches);
				if nargout == 2
					nmatches = length(statIndex);
					for n = 1:nmatches
						statObj(n) = obj.Statistics(statIndex(n));
					end
				end
			else
				statIndex = [];
				statObj = [];
			end
			
			return
		end
		%------------------------------------------------------------------------
		
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------		
		function statnames = GetStatisticNameList(obj)
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		%	Input Args:
		% 	Output Args:
		%------------------------------------------------------------------------
			if obj.Nstatistics == 0
				statnames = {};
				return
			end
		
			statnames = cell(obj.Nstatistics, 1);
			for s = 1:obj.Nstatistics
				statnames{s} = obj.Statistics(s).Name;
			end
			return
		end
		%------------------------------------------------------------------------
		
	end		%% END OF METHODS
end		%% END OF CLASS DEFINITION
