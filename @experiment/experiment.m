%---------------------------------------------------------------------
% experiment
%---------------------------------------------------------------------
% IntraHEKA Toolbox
%---------------------------------------------------------------------
% Class Definition
%---------------------------------------------------------------------
% 	Properties
%  		Sweeps						array of sweep objects holding sweep data
% 		Nsweeps						# of sweeps
% 		PGF							PGF object
% 		PUL							PUL object
% 		BaseName						base filename
% 		BasePath						base path
% 		DATfilename					.DAT filename
% 		PGFfilename					.PGF filename
% 		PULfilename					.PUL filename
% 		EXPfilename					.EXP filename
% 		Info							Experiment information
% 		PlotOpts						Plotting options
% 		isInitialized				0 if data are not loaded, 1 if they are
% 		Statistics					Statistic objects
% 		Nstatistics					# of statistic obejects
% 		
% 	Methods
% 		AddStatistic
% 		GetTrace
% 		BuildSweeps
% 		GetTraceLengthForCondition
% 		CheckInitAndCondition
% 		Initialize
% 		GetRawTrace
% 		MeanTraceForCondition
% 		GetResampledTrace
% 		PlotIndividualTracesForCondition
% 		GetSampleRate
% 		PlotTracesForCondition
% 		GetStatisticByName
% 		SetInfoFromLog
% 		GetStatisticNameList
% 		SpiketimesForCondition            
% 		GetStimParamForCondition
% 		GetStimulusForCondition
% 		GetSweepListForCondition
% 
%---------------------------------------------------------------------
% See also: sweep (class), readPGF, readPUL, readHEKA, readDAT
%---------------------------------------------------------------------

%---------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neomed.edu
%---------------------------------------------------------------------
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
%  3 May 2017 (SJS): 
%	 - added documentation
%	 - GetResampledTrace will now take a new sampling rate instead of 
% 		Decimation factor
%	 - original GetResampledTrace renames GetDecimatedTrace 
%---------------------------------------------------------------------
% TO DO:
%---------------------------------------------------------------------

%--------------------------------------------------------------------
% experiment class definition
%--------------------------------------------------------------------
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
		%--------------------------------------------------------------------
		%--------------------------------------------------------------------
		% Constructor
		%--------------------------------------------------------------------
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
			[base_path, base_name, ext] = fileparts(fullfile(base_path, ...
																base_name)); %#ok<ASGLU>
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
						warning('%s: unknown command %s', mfilename, ...
											upper(cmdstr));
				end
			end
		end
		%--------------------------------------------------------------------

		%--------------------------------------------------------------------
		%--------------------------------------------------------------------
		function Initialize(obj)
		%--------------------------------------------------------------------
		% initialize object (read in data)
		% must be called before doing most of the useful things in the 
		% experiment object.  
		%--------------------------------------------------------------------
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
		%--------------------------------------------------------------------
		
		%--------------------------------------------------------------------
		%--------------------------------------------------------------------
		function BuildSweeps(obj)
		%--------------------------------------------------------------------
		% builds the Sweeps object array property.  
		% This is usually called by the experiment.Initialize() method, but
		% can be called independently if you know what you're doing...
		%--------------------------------------------------------------------
			% use buildDSC function to get sweep information
			dsc = buildDSC(fullfile(obj.BasePath, obj.DATfilename), ...
															obj.PGF, obj.PUL.rr);
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
		%--------------------------------------------------------------------
		
		%--------------------------------------------------------------------
		%--------------------------------------------------------------------
		function [trace1, trace2] = GetRawTrace(obj, sweeplist)
		%--------------------------------------------------------------------
		% get raw (unscaled) trace(s)
		% returns [trace1, trace2], where trace1 is usually the stimulus, 
		% and trace2 is the electrode data
		%--------------------------------------------------------------------
			if ~obj.isInitialized
				warning('%s: object not initialized', mfilename);
				return
			end
			if any(~between(sweeplist, 1, obj.Nsweeps))
				error('%s: sweeplist out of bounds (Nsweeps = %d)', ...
																		mfilename, obj.Nsweeps);
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
				[trace1{n}, trace2{n}] = ...
								obj.Sweeps(sweeplist(n)).ReadData(datfile);
				trace1{n} = trace1{n};
				trace2{n} = trace2{n};
			end
		end
		%--------------------------------------------------------------------
		
		%--------------------------------------------------------------------
		%--------------------------------------------------------------------
		function [trace1, trace2] = GetTrace(obj, sweeplist)
		%--------------------------------------------------------------------
		% get scaled trace(s)
		% returns [trace1, trace2], where trace1 is usually the stimulus, 
		% and trace2 is the electrode data
		%--------------------------------------------------------------------
			if ~obj.isInitialized
				warning('%s: object not initialized', mfilename);
				return
			end
			
			if any(~between(sweeplist, 1, obj.Nsweeps))
				error('%s: sweeplist out of bounds (Nsweeps = %d)', ...
									mfilename, obj.Nsweeps);
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
				[trace1{n}, trace2{n}] = ...
											obj.Sweeps(sweeplist(n)).ReadData(datfile);
				trace1{n} = obj.Info.Scale(1) .* trace1{n} ./ obj.Info.Gain(1);
				trace2{n} = obj.Info.Scale(2) .* trace2{n} ./ obj.Info.Gain(2);
			end
		end
		%--------------------------------------------------------------------
		
		%--------------------------------------------------------------------
		%--------------------------------------------------------------------
		function [rate, trace1, trace2] = GetDecimatedTrace(obj, ...
																	sweeplist, decifactor)
		%--------------------------------------------------------------------
		% get scaled and decimated trace(s)
		% returns [rate, trace1, trace2], where trace1 is usually the stimulus, 
		% and trace2 is the electrode data. rate is new sampling rate
		%--------------------------------------------------------------------
			% check decimation factor
			if decifactor < 1
				error('%s: decifactor must be integer >= 1 (decifactor = %d)', ...
													mfilename, decifactor);
			end
			% get the trace(s) in raw format
			[trace1, trace2] = obj.GetTrace(sweeplist);
			% if decifactor == 1, nothing to do!
			if decifactor == 1
				% Rate from Sweeps is actually dt
				rate = 1./obj.Sweeps(sweeplist(1)).Rate;
				return
			end
			% # of sweeps
			ns = length(sweeplist);
			% determine original sample rate - this is a bit of kludge, under
			% the possibility that sample rate might vary across sweeps, 
			% which is unlikely...  nevertheless...
			tmprate = zeros(size(sweeplist));
			for n = 1:ns
				tmprate(n) = 1./obj.Sweeps(sweeplist(n)).Rate;
			end
			% use the lowest value of the found rates (actually, sample
			% intervals, in seconds) and multiply by decifactor to get the new
			% rate
			rate = min(tmprate) / decifactor;
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
		%--------------------------------------------------------------------

		%--------------------------------------------------------------------
		%--------------------------------------------------------------------
		function [trace1, trace2, varargout] = GetResampledTrace(obj, ...
																	sweeplist, newFs)
		%--------------------------------------------------------------------
		% get scaled and decimated trace(s)
		% returns [trace1, trace2, oldFs], where trace1 is usually the stimulus, 
		% and trace2 is the electrode data. oldFs is original sampling rate
		% uses 	y = resample(x, tx, newFs) method
		%--------------------------------------------------------------------
			% check newFs
			if newFs < 1
				error('%s: newFs must be >= 1', mfilename, newFs);
			end
			% get the trace(s) in raw format
			[trace1, trace2] = obj.GetTrace(sweeplist);
			% # of sweeps
			nSweeps = length(sweeplist);
			% determine original sample rate - this is a bit of kludge, under
			% the possibility that sample rate might vary across sweeps, 
			% which is unlikely...  nevertheless...
			oldFs = zeros(size(sweeplist));
			nsamples = zeros(size(sweeplist));
			for n = 1:nSweeps
				% Rate from Sweeps is actually dt
				oldFs(n) = 1/obj.Sweeps(sweeplist(n)).Rate;
				nsamples(n) = obj.Sweeps(sweeplist(n)).Nsamples;
			end
			% do a quick check on rates
			if any(oldFs(1) ~= oldFs)
				warning('some sweeps have different sampling rates!!!!')
				RATEMATCH = 0;
			else
				RATEMATCH = 1;
				oldFs = oldFs(1);
			end
			% and check samples
			if any(nsamples(1) ~= nsamples)
				warning('some sweeps have different lengths!!!!!');
				NSAMPMATCH = 0;
			else
				NSAMPMATCH = 1;
				nsamples = nsamples(1);
			end
			% if all rates match and newFs is same, we're done
			if RATEMATCH && (oldFs(1) == newFs)
				if nargout > 2
					varargout{1} = oldFs;
				end
				return
			end
			
			% now, resample, taking into account the possible mismatches
			if RATEMATCH && NSAMPMATCH
				% this is easiest, and fastest, situation, since all sweeps are
				% the same sample rate and length
				% only need to create 1 time vector
				% build time base for orig data
				t_orig = (0:(nsamples(1) - 1)) * (1/oldFs(1))';
				% resample original data
				for n = 1:nSweeps
					trace1{n} = resample(trace1{n}, t_orig, newFs, 'pchip');
					trace2{n} = resample(trace2{n}, t_orig, newFs, 'pchip');
				end
			else
				for n = 1:nSweeps
					% build time base for orig data
					t_orig = (0:(nsamples(n) - 1)) * (1/oldFs(n))';
					% resample original data
					trace1{n} = resample(trace1{n}, t_orig, newFs, 'pchip');
					trace2{n} = resample(trace2{n}, t_orig, newFs, 'pchip');
				end				
			end
		end
		%--------------------------------------------------------------------
		
		
		%--------------------------------------------------------------------
		%--------------------------------------------------------------------
		function [Fs, varargout] = GetSampleRate(obj)
			% determine original sample rate - this is a bit of kludge, under
			% the possibility that sample rate might vary across sweeps, 
			% which is unlikely...  nevertheless...
			tmprate = zeros(obj.Nsweeps, 1);
			for n = 1:obj.Nsweeps
				tmprate(n) = obj.Sweeps(n).Rate;
			end
			% use the lowest value of the found rates (actually, sample
			% intervals, in seconds)
			Fs = 1/min(tmprate);
			if nargout > 1
				varargout{1} = unique(tmprate.^-1);
			end
		end
		%--------------------------------------------------------------------
		%--------------------------------------------------------------------
		
		%--------------------------------------------------------------------
		%--------------------------------------------------------------------
		function SetInfoFromLog(obj, logData)
		%--------------------------------------------------------------------
		% set stimulus information
		% apply log info to sweeps
		%--------------------------------------------------------------------
			% make sure object is initialized
			if ~obj.isInitialized
				warning('%s: object not initialized', mfilename);
				return
			end
			% first, need to go through and combine the year, month, day unit
			% fields
			nLog = length(logData);
			unitNames = cell(nLog, 1);
			for l = 1:nLog
				L = logData(l);
				unitNames{l} = sprintf('%s-%s-%s-%s', L.year, L.month, ...
													L.day, L.unit_number);
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
					Stimulus(n).Ntrials = 1 + Stimulus(n).TrialEnd - ...
													Stimulus(n).TrialStart;
					Stimulus(n).AuditoryStimulus = L.auditory_stim;
					Stimulus(n).OtherStimulus = L.other_stim;
					Stimulus(n).Comments = L.comments;
				end
				% store Stimulus in object
				obj.Info.Stimulus = Stimulus;
			end
			
		end
		%--------------------------------------------------------------------
		
		%--------------------------------------------------------------------
		%--------------------------------------------------------------------
		function [t2avg, t2std, varargout] = MeanTraceForCondition(obj, ...
																		Condition, varargin)
		%--------------------------------------------------------------------
		% returns average trace for specific Condition
		%--------------------------------------------------------------------
		% [traceavg, tracestd, stimtrace, samprate, traces] = 
		%								MeanTraceForCondition(Condition,	
		% 															 'NoSpikes',
		%															 'Resample', <new rate>,
		%																 OR
		% 															 'Decimate', <deci factor>)
		%--------------------------------------------------------------------
			% make sure object is initialized
			if ~obj.CheckInitAndCondition(Condition)
				error('%s: Condition %d not found or class is not initialized', ...
								mfilename, Condition);
			end
			% check input args
			DESPIKE = false;
			RESAMPLE = false;
			DECIMATE = false;
			optargin = size(varargin,2);
% 			stdargin = nargin - optargin;
			if optargin
				% parse varargin args
				n = 1;
				while n < optargin
					if ischar(varargin{n})
						if strcmpi(varargin{n}, 'NoSpikes')
							DESPIKE = true;
							n = n + 1;
						elseif strcmpi(varargin{n}, 'Resample')
							RESAMPLE = true;
							new_rate = varargin{n+1};
							n = n+2;
						elseif strcmpi(varargin{n}, 'Decimate')
							DECIMATE = true;
							deci_factor = varargin{n+1};
							n = n+1;
						else
							error('%s: unknown option %s', mfilename, varargin{n});
						end
					end
				end
			end
			if RESAMPLE && DECIMATE
				warning('%s: cannot both resample and decimate!', mfilename);
				RESAMPLE = false;
				DECIMATE = false;
			end
			% get the sweeps for the indicated Condition (resample or decimate
			% as needed)
			obj.GetStimulusForCondition(Condition);
			sweeplist = obj.GetSweepListForCondition(Condition);
			if RESAMPLE
				[trace1, trace2] = obj.GetResampledTrace(sweeplist, new_rate);
				Fs = new_rate;
			elseif DECIMATE
				[Fs, trace1, trace2] = obj.GetDecimatedTrace(sweeplist, ...
																					deci_factor);
			else
				[trace1, trace2] = obj.GetTrace(sweeplist);
				Fs = obj.GetSampleRate;
			end
			% store first trace1 data, then clear trace1 to save space
			t1 = trace1{1};
			clear trace1
			% despike traces if needed
			if DESPIKE
				for s = 1:length(sweeplist)
					trace2{s} = deSpike(trace2{s}, Fs);
				end
			end
			% compute mean of trace2 cell array
			t2avg = mean(cell2mat(trace2'), 2);
			t2std = std(cell2mat(trace2'), 0, 2);
			% assign outputs
			if nargout >= 3
				varargout{1} = t1;
			end
			if nargout >= 4
				varargout{2} = Fs;
			end
			if nargout == 5
				varargout{3} = trace2;
			end			
		end
		%--------------------------------------------------------------------
		
		
		%--------------------------------------------------------------------
		%--------------------------------------------------------------------
		function [spiket, nspikes] = SpiketimesForCondition(obj, ...
																		Condition, varargin)
		%--------------------------------------------------------------------
		% returns spike times for specific Condition
		%--------------------------------------------------------------------
			%-------------------------------------
			% default settings
			%-------------------------------------
			% threshold for spikes, V
			spikethresh = -0.015;
			% spike refractory period, ms
			refractorytime = 2;			
			%-------------------------------------
			% make sure object is initialized
			%-------------------------------------
			if ~obj.CheckInitAndCondition(Condition)
				error('%s: Condition %d not found or class is not initialized', ...
								mfilename, Condition);
			end
			%-------------------------------------
			% check inputs
			%-------------------------------------
			if ~isempty(varargin)
				optargs = length(varargin);
				% parse varargin args
				n = 1;
				while n < optargs
					switch upper(varargin{n})
						case 'SPIKETHRESHOLD'
							spikethresh = varargin{n+1};
							n = n + 2;
						case 'REFRACTORYTIME'
							refractorytime = varargin{n+1};
							n = n + 2;
						otherwise
							error('%s: bad arg %s', mfilename, varargin{n});
					end
				end
			end
			%-------------------------------------
			% get the sweeps for the indicated Condition
			%-------------------------------------
			obj.GetStimulusForCondition(Condition);		
			sweeplist = obj.GetSweepListForCondition(Condition);
			[~, traces] = obj.GetTrace(sweeplist);
			%-------------------------------------			
			% get sample rate
			%-------------------------------------
			Fs = obj.GetSampleRate;
			%-------------------------------------
			% find spike times (convert bins to seconds)
			%-------------------------------------			
			spiket = cell(size(traces));
			nspikes = zeros(size(traces));
			for s = 1:length(sweeplist)
				spiket{s} = spikeschmitt2(force_row(traces{s}) - spikethresh, ...
													0, refractorytime, Fs);
				spiket{s} = spiket{s} ./ Fs;
				nspikes(s) = length(spiket{s});
			end
		end
		%--------------------------------------------------------------------
		
		%--------------------------------------------------------------------
		%--------------------------------------------------------------------
		function tlen = GetTraceLengthForCondition(obj, Condition) 
		%--------------------------------------------------------------------
		% Returns trace lengths for a given condition
		%--------------------------------------------------------------------
			% make sure object is initialized
			if ~obj.CheckInitAndCondition(Condition)
				tlen = [];
				return
			end
			% get the sweeps for the indicated Condition
			sinfo = obj.GetStimulusForCondition(Condition);		 %#ok<NASGU>
			sweeplist = obj.GetSweepListForCondition(Condition);
			tlen = zeros(length(sweeplist), 1);
			% get # of samples for each sweep
			for t = 1:length(tlen)
				tlen(t) = obj.Sweeps(sweeplist(t)).Nsamples;
			end
		end
		%--------------------------------------------------------------------
		
		%--------------------------------------------------------------------
		%--------------------------------------------------------------------
		function sweeplist = GetSweepListForCondition(obj, Condition) 
		%--------------------------------------------------------------------
		% returns sweep list (indices) for given Condition
		%--------------------------------------------------------------------
			% make sure object is initialized
			if ~obj.CheckInitAndCondition(Condition)
				sweeplist = [];
				return
			end
			% get the sweeps for the indicated Condition
			sinfo = obj.GetStimulusForCondition(Condition);		
			sweeplist = sinfo.TrialStart : sinfo.TrialEnd;
		end
		%--------------------------------------------------------------------
		
		%--------------------------------------------------------------------
		%--------------------------------------------------------------------
		function s = GetStimulusForCondition(obj, Condition)
		%--------------------------------------------------------------------
		% Method to get Stimulus information for specific condition
		%--------------------------------------------------------------------
			% make sure object is initialized
			if ~obj.CheckInitAndCondition(Condition)
				s = [];
				return
			else
				s = obj.Info.Stimulus(Condition);
			end
		end
		%--------------------------------------------------------------------

		%--------------------------------------------------------------------
		%--------------------------------------------------------------------
		function out = CheckInitAndCondition(obj, Condition)
		%--------------------------------------------------------------------
		% internal method to check initialization status and validity of 
		% Condition value
		%--------------------------------------------------------------------
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
		%--------------------------------------------------------------------
		
		%--------------------------------------------------------------------
		%--------------------------------------------------------------------		
		function stim = GetStimParamForCondition(obj, Condition, varargin) 
		%--------------------------------------------------------------------
		% Returns stimulus information for a given condition
		% Note that these are parameters that are computed directly from the 
		% stimulus trace and may be inaccurate!!!
		%--------------------------------------------------------------------
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
		%--------------------------------------------------------------------
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
			obj.GetStimulusForCondition(Condition);
			sweeplist = obj.GetSweepListForCondition(Condition);
			[rate, trace1, ~] = obj.GetDecimatedTrace(sweeplist(1), DFACT);
			% then, compute envelope using Hilbert transform and find start
			% and end
			env = abs(hilbert(1000*trace1));
			stim.start = 1000 * rate * find(env > THRESHOLD, 1, 'first');
			stim.end = 1000 * rate * find(env > THRESHOLD, 1, 'last');
		end
		%--------------------------------------------------------------------
		
		
		%--------------------------------------------------------------------
		%--------------------------------------------------------------------
		function [statIndex, statObj] = AddStatistic(obj, Statname, Statvals)
		%--------------------------------------------------------------------
		%--------------------------------------------------------------------
		%	Input Args:
		% 	Output Args:
		%--------------------------------------------------------------------
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
		%--------------------------------------------------------------------
		
		%--------------------------------------------------------------------
		%--------------------------------------------------------------------		
		function [statIndex, statObj] = GetStatisticByName(obj, Statname, ...
																				CompareMethod)
		%--------------------------------------------------------------------
		%--------------------------------------------------------------------
		%	Input Args:
		% 	Output Args:
		%--------------------------------------------------------------------
			% check inputs
			if obj.Nstatistics == 0
				statIndex = 0;
				statObj = [];
				return
			end			
			if ~exist('CompareMethod', 'var')
				CompareMethod = 'STRCMPI';
			end
			% add stat
			statnamelist = obj.GetStatisticNameList;
			switch upper(CompareMethod)
				case 'STRCMP'
					namematches = strcmp(Statname, statnamelist);
				case 'STRCMPI'
					namematches = strcmpi(Statname, statnamelist);
				case 'STRNCMP'
					namematches = strncmp(Statname, statnamelist, ...
														length(Statname));
				case 'STRNCMPI'
					namematches = strncmpi(Statname, statnamelist, ...
														length(Statname));
				otherwise
					error('%s: unknown stat string compare method %s', ...
																mfilename, CompareMethod);
			end
			% check for matches
			if any(namematches)
				statIndex = find(namematches);
				if nargout == 2
					nmatches = length(statIndex);
					for n = 1:nmatches
						statObj(n) = obj.Statistics(statIndex(n)); %#ok<AGROW>
					end
				end
			else
				statIndex = [];
				statObj = [];
			end
			% done
			return
		end
		%--------------------------------------------------------------------
		
		%--------------------------------------------------------------------
		%--------------------------------------------------------------------		
		function statnames = GetStatisticNameList(obj)
		%--------------------------------------------------------------------
		%--------------------------------------------------------------------
		%	Input Args:
		% 	Output Args:
		%--------------------------------------------------------------------
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
		%--------------------------------------------------------------------
		
	end		%% END OF METHODS
end		%% END OF CLASS DEFINITION
