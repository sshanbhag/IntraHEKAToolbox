function varargout = HEKAview(varargin)
% HEKAVIEW MATLAB code for HEKAview.fig
%      HEKAVIEW, by itself, creates a new HEKAVIEW or raises the existing
%      singleton*.
%
%      H = HEKAVIEW returns the handle to a new HEKAVIEW or the handle to
%      the existing singleton*.
%
%      HEKAVIEW('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in HEKAVIEW.M with the given input arguments.
%
%      HEKAVIEW('Property','Value',...) creates a new HEKAVIEW or raises
%      the existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before HEKAview_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to HEKAview_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help HEKAview

% Last Modified by GUIDE v2.5 08-May-2017 14:18:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @HEKAview_OpeningFcn, ...
                   'gui_OutputFcn',  @HEKAview_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT
%-------------------------------------------------------------------------- 

%*****************************************************************************
%*****************************************************************************
%*****************************************************************************
%	GUI functions/callbacks
%*****************************************************************************
%*****************************************************************************
%*****************************************************************************

%-------------------------------------------------------------------------- 
%-------------------------------------------------------------------------- 
% --- Executes just before HEKAview is made visible.
%-------------------------------------------------------------------------- 
function HEKAview_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>
% This function has no output args, see OutputFcn.
%	hObject    handle to figure
%	eventdata  reserved - to be defined in a future version of MATLAB
%	handles    structure with handles and user data (see GUIDATA)
%	varargin   command line arguments to HEKAview (see VARARGIN)

	% Choose default command line output for HEKAview
	handles.output = hObject;
	% Update handles structure
	guidata(hObject, handles);
	initialize_gui(hObject, handles, false);
	% UIWAIT makes HEKAview wait for user response (see UIRESUME)
	% uiwait(handles.figure1);
%-------------------------------------------------------------------------- 
%-------------------------------------------------------------------------- 

%-------------------------------------------------------------------------- 
%-------------------------------------------------------------------------- 
% --- Outputs from this function are returned to the command line.
%-------------------------------------------------------------------------- 
function varargout = HEKAview_OutputFcn(hObject, eventdata, handles)
% 	varargout  cell array for returning output args (see VARARGOUT);
% 	hObject    handle to figure
% 	eventdata  reserved - to be defined in a future version of MATLAB
% 	handles    structure with handles and user data (see GUIDATA)

	% Get default command line output from handles structure
	varargout{1} = handles.output;
%-------------------------------------------------------------------------- 
%-------------------------------------------------------------------------- 


%**************************************************************************
%**************************************************************************
%**************************************************************************
% CONDITION and SWEEP panel callbacks
%**************************************************************************
%**************************************************************************
%**************************************************************************
%-------------------------------------------------------------------------- 
%-------------------------------------------------------------------------- 
function Condition_ctrl_Callback(hObject, eventdata, handles)
	%----------------------------------------------------------
	% make sure there's something there
	%----------------------------------------------------------
	if isempty(read_ui_str(hObject))
		return
	elseif ~handles.E.isInitialized
		warning('Experiment Object is not initialized!')
		fprintf('\n¡¡Please load data!!\n\n');
		return
	end
	handles.Values.CurrentCondition = read_ui_val(hObject);
	%----------------------------------------------------------
	% update Sweep list box
	%----------------------------------------------------------
	% first get the list of sweeps for current condition (and # of sweeps)
	sweeplist = handles.E.GetSweepListForCondition( ...
												handles.Values.CurrentCondition);
	handles.Values.Nsweeps = length(sweeplist);
	handles.Values.Sweeplist = sweeplist;
	% set the current sweep to -1 (all sweeps)
	handles.Values.CurrentSweep = -1;
	% create text list of sweeps
	cstr = cell(handles.Values.Nsweeps + 1, 1);
	for n = 1:handles.Values.Nsweeps
		cstr{n} = sprintf('%d', sweeplist(n));
	end
	cstr{end} = 'All';
	set(handles.Sweep_ctrl, 'String', cstr);		
	% update Sweep control to current sweep (all sweeps)
	update_ui_val(handles.Sweep_ctrl, handles.Values.Nsweeps + 1);
	% get # of points for this condition
	npts = min(handles.E.GetTraceLengthForCondition(...
												handles.Values.CurrentCondition));
	% get the stimulus info
	sinfo = handles.E.Info.Stimulus(handles.Values.CurrentCondition);
	% get the information for the first sweep
	firstsweep = sinfo.TrialStart;
	% use sample rate from this sweep
	rate = handles.E.Sweeps(firstsweep).Rate;
	% set Tmin and Tmax values
	handles.Values.Tmin = 0;
	handles.Values.Tmax = npts * rate * 1000;
	guidata(hObject, handles);
	% update gui elements
	update_ui_str(handles.Tmin_ctrl, handles.Values.Tmin); 
	update_ui_str(handles.Tmax_ctrl, handles.Values.Tmax); 
	update_gui(hObject, handles);
%-------------------------------------------------------------------------- 
%-------------------------------------------------------------------------- 
function Sweep_ctrl_Callback(hObject, eventdata, handles)
	index = read_ui_val(hObject);
	if index <= handles.Values.Nsweeps
		handles.Values.CurrentSweep = handles.Values.Sweeplist(index);
	else
		% -1 corresponds to ALL sweeps
		handles.Values.CurrentSweep = -1;
	end
	guidata(hObject, handles);
	update_gui(hObject, handles);
%-------------------------------------------------------------------------- 
%-------------------------------------------------------------------------- 
%**************************************************************************
%**************************************************************************
%**************************************************************************


%**************************************************************************
%**************************************************************************
%**************************************************************************
%-------------------------------------------------------------------------- 
%--------------------------------------------------------------------------
% LOAD DATA
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function buttonLoadData_Callback(hObject, eventdata, handles)
	LoadData(hObject, eventdata, handles);
%-------------------------------------------------------------------------- 
%**************************************************************************
%**************************************************************************
%**************************************************************************


%**************************************************************************
%**************************************************************************
%**************************************************************************
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% CALCULATE PANEL CALLBACKS
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function buttonMean_Callback(hObject, eventdata, handles)
%--------------------------------------------------------------------------
% Calculate mean trace
%--------------------------------------------------------------------------
	%-----------------------------------------------
	% despike?
	%-----------------------------------------------
	if handles.Values.RemoveSpikesFromMean
		disp('Removing Spikes for mean calculation');
		[sweep_mean, sweep_std, stim, Fs, sweeps] = ...
										handles.E.MeanTraceForCondition( ...
												handles.Values.CurrentCondition, ...
												'NoSpikes', ...
												'Decimate', handles.Values.Decifactor);
	else
		[sweep_mean, sweep_std, stim, Fs, sweeps] = ...
										handles.E.MeanTraceForCondition( ...
												handles.Values.CurrentCondition, ...
												'Decimate', handles.Values.Decifactor);
	end
	%-----------------------------------------------
	% plot mean and error bars
	%-----------------------------------------------
	% time vector
	tvec = 1000 * ((1:length(sweep_mean)) - 1) * (1/Fs);
	% create new fig
	figure
	% plot mean and area (convert mean and std to mV from V
	shadedErrorBar(tvec, 1000*sweep_mean, 1000.*sweep_std, ...
									{'MarkerFaceColor', [0 0.4470 0.7410]})
	xlabel('Time (ms)')
	ylabel('mV');
	fname = sprintf('%s_%d-%d.mat', ...
								handles.Files.basename, ...
								handles.Values.Sweeplist(1), ...
								handles.Values.Sweeplist(end) );
	title(fname, 'Interpreter', 'none');
	% plot stimulus at top of plot
	yminmax = ylim;
	stim2 = yminmax(2) - 10*( diff(yminmax)/100)*normalize(stim);
	hold on
		plot(tvec, stim2, 'Color', 0.5 * [1 1 1])
	hold off
	ylim([yminmax(1) max(stim2)]);
	if handles.Values.OverlaySweeps
		hold on
		for n = 1:length(sweeps)
			plot(tvec, 1000*sweeps{n}, 'Color', 0.25 * [1 1 1]);
		end
	end	
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function OverlaySweeps_ctrl_Callback(hObject, eventdata, handles)
	% Controls whether individual spikes will be plotted along with computed
	% mean trace
	handles.Values.OverlaySweeps = read_ui_val(hObject);
	guidata(hObject, handles);
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function RemoveSpikesFromMean_ctrl_Callback(hObject, eventdata, handles)
	% controls whether spikes will included or removed from mean calculation
	handles.Values.RemoveSpikesFromMean = read_ui_val(hObject);
	guidata(hObject, handles);
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% --- Executes on button press in buttonSpikeTimes.
%--------------------------------------------------------------------------
function buttonSpikeTimes_Callback(hObject, eventdata, handles)
	%-----------------------------------------------
	% get threshold, holdoff values
	%-----------------------------------------------
	% threshold for spikes, V
	spikethresh = 0.001 * handles.Values.SpikeThreshold;
	% spike refractory period, ms
	refractorytime = handles.Values.HoldoffTime;
	%-----------------------------------------------
	% check init
	%-----------------------------------------------
	if ~handles.E.isInitialized
		warning('%s: Not initialized!', mfilename);
		return
	end
	%-----------------------------------------------
	% get sweeps
	%-----------------------------------------------
	[spiket, nspikes] = ...
					handles.E.SpiketimesForCondition( ...
								handles.Values.CurrentCondition, ...
								'SPIKETHRESHOLD', spikethresh, ...
								'REFRACTORYTIME', refractorytime ); %#ok<ASGLU>
			
	%-----------------------------------------------
	% plot rasters
	%-----------------------------------------------
	figure
	rasterplot(spiket);
% 	%-----------------------------------------------
% 	% get output filename from user
% 	%-----------------------------------------------
% 	% build default file name
% 	fname = sprintf('%s_%d-%d_spiketimes.mat', ...
% 								fullfile(handles.Files.rawpath, ...
% 												handles.Files.basename), ...
% 								handles.Values.Sweeplist(1), ...
% 								handles.Values.Sweeplist(end) );
% 	[filename, pathname] = uiputfile('*.mat', 'Export As', fname);
% 	% if user cancelled, abort
% 	if isequal(filename, 0) || isequal(pathname, 0)
% 		disp('Cancelled Spiketime File..')
% 		return
% 	else
% 		save(fullfile(pathname, filename), 'spiket', 'nspikes', '-MAT');
% 	end
	guidata(hObject, handles);
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function SpikeThreshold_ctrl_Callback(hObject, eventdata, handles)
%---------------------------------------
% set Spike Threshold
%---------------------------------------
	newVal = read_ui_str(handles.SpikeThreshold_ctrl, 'n');
	if ~isnumeric(newVal) || isempty(newVal)
		errordlg('Spike Threshold must be a number', ...
					'HEKAview: Threshold error');
		update_ui_str(handles.SpikeThreshold_ctrl, handles.Values.SpikeThreshold);
	elseif ~between(newVal, -1000, 1000)
		errordlg('Spike Threshold must be between -1000 and 1000 mV', ...
					'HEKAview: Threshold error');
		update_ui_str(handles.SpikeThreshold_ctrl, handles.Values.SpikeThreshold);
	else
		handles.Values.SpikeThreshold = newVal;
	end
	guidata(hObject, handles);
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function SpikeHoldoff_ctrl_Callback(hObject, eventdata, handles)
%---------------------------------------
% set Spike Holdoff
%---------------------------------------
	newVal = read_ui_str(handles.SpikeHoldoff_ctrl, 'n');
	if ~isnumeric(newVal) || isempty(newVal)
		errordlg('Spike Holdoff value must be a number', ...
					'HEKAview: Holdoff error');
		update_ui_str(handles.SpikeHoldoff_ctrl, handles.Values.HoldoffTime);
	elseif ~between(newVal, 0, 1000)
		errordlg('Spike holdoff time must be between 0 and 1000 ms', ...
					'HEKAview: Holdoff error');
		update_ui_str(handles.SpikeHoldoff_ctrl, handles.Values.HoldoffTime);
	else
		handles.Values.HoldoffTime = newVal;
	end
	guidata(hObject, handles);
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%**************************************************************************
%**************************************************************************
%**************************************************************************


%**************************************************************************
%**************************************************************************
%**************************************************************************
%--------------------------------------------------------------------------
% LOG FILE PANEL CALLBACKS
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%-------------------------------------------------------------------------- 
function buttonLoadLogFile_Callback(hObject, eventdata, handles)
	%-----------------------------------------------
	% get file from user
	% 	logpath = '/Users/sshanbhag/Work/Data/Mouse/IntracellularAmygdala/RawData';
	% 	logfile = 'IntracellDataLog.csv';
	% 	logname = fullfile(handles.Files.logpath, handles.Files.logfile);
	%-----------------------------------------------
	[filename, pathname] = uigetfile('*.csv', ...
						'Select experiment log file', handles.Files.logpath);
	% if user cancelled, abort
	if isequal(filename, 0) || isequal(pathname, 0)
		disp('Cancelled Load File..')
		return
	end	
	handles.Files.logpath = pathname;
	handles.Files.logfile = filename;
	handles.Files.logname = fullfile(handles.Files.logpath, ...
													handles.Files.logfile);
	set(handles.textLogFile, 'String', handles.Files.logname);
	guidata(hObject, handles);
%-------------------------------------------------------------------------- 
%-------------------------------------------------------------------------- 
%**************************************************************************
%**************************************************************************
%**************************************************************************


%**************************************************************************
%**************************************************************************
%**************************************************************************
% --- Export Callbacks
%**************************************************************************
%**************************************************************************
%**************************************************************************
%-------------------------------------------------------------------------- 
%-------------------------------------------------------------------------- 
%-------------------------------------------------------------------------- 
function ExportSweeps_ctrl_Callback(hObject, eventdata, handles)
	%-----------------------------------------------
	% check init
	%-----------------------------------------------
	if ~handles.E.isInitialized
		warning('%s: Not initialized!', mfilename);
		return
	end
	%-----------------------------------------------
	% export mode
	%-----------------------------------------------
	fprintf('Exporting data as %s\n', handles.Values.Export.ExportFormat);
	%-----------------------------------------------
	% downsample? get sweeps, info
	%-----------------------------------------------
	Fs = handles.E.GetSampleRate;
	newFs = handles.Values.Export.SampleRate;
	if handles.Values.Export.ResampleData && (Fs ~= newFs)
		[Stimuli, Sweeps] = handles.E.GetResampledTrace(...
													handles.Values.Sweeplist, ...
													newFs);
		fprintf('Original SampleRate: %d\n', Fs);
		fprintf('Export SampleRate: %d\n', newFs);
		Fs = newFs;
	else
		[Stimuli, Sweeps] = handles.E.GetTrace(handles.Values.Sweeplist);		
	end
	nStimuli = length(Stimuli);
	nSweeps = length(Sweeps);
	%-----------------------------------------------
	% despike?
	%-----------------------------------------------
	if handles.Values.Export.RemoveSpikes
		for s = 1:length(Sweeps)
			[Sweeps{s}, ~] = deSpike(Sweeps{s}, Fs);
		end
	end
	%-----------------------------------------------
	% stimulus info
	%-----------------------------------------------
	AuditoryStimulus = handles.AuditoryStim_text.String;
	OtherStimulus = handles.OtherStim_text.String;
	Comments = handles.Comments_text.String;
	
	%-----------------------------------------------
	% Export!
	%-----------------------------------------------
	if strcmpi(handles.Values.Export.ExportFormat, 'CSV')
		%-----------------------------------------------
		% export as csv
		%-----------------------------------------------
		%-----------------------------------------------
		% get output filename from user
		%-----------------------------------------------
		% build default file name
		fname = sprintf('%s_%d-%d.csv', ...
									fullfile(handles.Files.rawpath, ...
													handles.Files.basename), ...
									handles.Values.Sweeplist(1), ...
									handles.Values.Sweeplist(end) );
		[filename, pathname] = uiputfile('*.csv', 'Export As Text File', fname);
		% if user cancelled, abort
		if isequal(filename, 0) || isequal(pathname, 0)
			disp('Cancelled Export File..')
			return
		end
		% check on sweeps
		if nStimuli ~= nSweeps
			error('%s: mismatch in sweeps and stimuli', mfilename);
		end
		%-----------------------------------------------
		% export to CSV
		%-----------------------------------------------
		fp = fopen(fullfile(pathname, filename), 'wt');
		fprintf(fp, 'AuditoryStim:,%s,\n', AuditoryStimulus);
		fprintf(fp, 'OtherStim:,%s,\n', OtherStimulus);
		fprintf(fp, 'Comments:,%s,\n', Comments);
		fprintf(fp, 'SampleRate:,%f,\n', Fs);
		fprintf(fp, 'Data:,\n');
		for n = 1:length(Sweeps{1})
			for s = 1:nSweeps
				fprintf(fp, '%f,%f,', Stimuli{s}(n), Sweeps{s}(n));
			end
			fprintf(fp, '\n');
		end
		fclose(fp);
	else
		%-----------------------------------------------
		% export as MAT
		%-----------------------------------------------
		%-----------------------------------------------
		% get output filename from user
		%-----------------------------------------------
		% build default file name
		fname = sprintf('%s_%d-%d.mat', ...
									fullfile(handles.Files.rawpath, ...
													handles.Files.basename), ...
									handles.Values.Sweeplist(1), ...
									handles.Values.Sweeplist(end) );
		[filename, pathname] = uiputfile('*.mat', 'Export As', fname);
		% if user cancelled, abort
		if isequal(filename, 0) || isequal(pathname, 0)
			disp('Cancelled Export File..')
			return
		end
		%-----------------------------------------------
		% export to MAT
		%-----------------------------------------------
		SweepList = handles.Values.SweepList; %#ok<NASGU>
		save(fullfile(pathname, filename), 'Stimuli', 'Sweeps', 'Fs', ...
									'AuditoryStimulus', 'OtherStimulus', 'Comments', ...
									'nStimuli', 'nSweeps', 'Sweeplist', '-MAT');
	end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function exportSpikeTimes_ctrl_Callback(hObject, eventdata, handles)
	%-----------------------------------------------
	% get threshold, holdoff values
	%-----------------------------------------------
	% threshold for spikes, V
	spikethresh = 0.001 * handles.Values.SpikeThreshold;
	% spike refractory period, ms
	refractorytime = handles.Values.HoldoffTime;
	%-----------------------------------------------
	% check init
	%-----------------------------------------------
	if ~handles.E.isInitialized
		warning('%s: Not initialized!', mfilename);
		return
	end
	%-----------------------------------------------
	% get sweeps
	%-----------------------------------------------
	[spiket, nspikes] = ...
					handles.E.SpiketimesForCondition( ...
								handles.Values.CurrentCondition, ...
								'SPIKETHRESHOLD', spikethresh, ...
								'REFRACTORYTIME', refractorytime );
	Fs = handles.E.GetSampleRate;
	%-----------------------------------------------
	% plot rasters
	%-----------------------------------------------
	figure
	rasterplot(spiket);
	%-----------------------------------------------
	% animal and stimulus info
	%-----------------------------------------------
	Animal = handles.E.Info.Animal;
	Depth = handles.E.Info.Depth;
	AuditoryStimulus = handles.AuditoryStim_text.String;
	OtherStimulus = handles.OtherStim_text.String;
	Comments = handles.Comments_text.String;
	%-----------------------------------------------
	% save
	%-----------------------------------------------
	if strcmpi(handles.Values.Export.ExportFormat, 'CSV')
		%-----------------------------------------------
		% export as csv
		%-----------------------------------------------
		%-----------------------------------------------
		% get output filename from user
		%-----------------------------------------------
		% build default file name
		fname = sprintf('%s_%d-%d_spiketimes.csv', ...
									fullfile(handles.Files.rawpath, ...
													handles.Files.basename), ...
									handles.Values.Sweeplist(1), ...
									handles.Values.Sweeplist(end) );
		[filename, pathname] = uiputfile('*.csv', 'Export As Text File', fname);
		% if user cancelled, abort
		if isequal(filename, 0) || isequal(pathname, 0)
			disp('Cancelled Spiketime File..')
			return
		end
		%-----------------------------------------------
		% export to CSV
		%-----------------------------------------------
		fp = fopen(fullfile(pathname, filename), 'wt');
		fprintf(fp, 'Animal:,%s,\n', Animal);
		fprintf(fp, 'Depth:,%d,\n', Depth);
		fprintf(fp, 'AuditoryStim:,%s,\n', AuditoryStimulus);
		fprintf(fp, 'OtherStim:,%s,\n', OtherStimulus);
		fprintf(fp, 'Comments:,%s,\n', Comments);
		fprintf(fp, 'SampleRate:,%f,\n', Fs);
		fprintf(fp, 'Sweep#,SpikeTimes:,\n');
		nSweeps = length(spiket);
		for s = 1:nSweeps
			fprintf(fp, '%d,', handles.Values.Sweeplist(s));
			for t = 1:nspikes(s)
				fprintf(fp, '%f,', spiket{s}(t));
			end
			fprintf(fp, '\n');
		end
		fclose(fp);		
	else
		%-----------------------------------------------
		% export as MAT
		%-----------------------------------------------
		%-----------------------------------------------
		% get output filename from user
		%-----------------------------------------------
		% build default file name
		fname = sprintf('%s_%d-%d_spiketimes.mat', ...
									fullfile(handles.Files.rawpath, ...
													handles.Files.basename), ...
									handles.Values.Sweeplist(1), ...
									handles.Values.Sweeplist(end) );
		[filename, pathname] = uiputfile('*.mat', 'Export As', fname);
		% if user cancelled, abort
		if isequal(filename, 0) || isequal(pathname, 0)
			disp('Cancelled Spiketime File..')
			return
		end
		SweepList = handles.Values.Sweeplist; %#ok<NASGU>
		save(fullfile(pathname, filename), 'spiket', 'nspikes', 'Fs', ...
									'AuditoryStimulus', 'OtherStimulus', 'Comments', ...
									'nStimuli', 'nSweeps', 'Sweeplist', '-MAT');
	end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function exportRemoveSpikes_ctrl_Callback(hObject, eventdata, handles)
	val = read_ui_val(hObject);
	handles.Values.Export.RemoveSpikes = val;
	guidata(hObject, handles);
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function exportSampleRate_ctrl_Callback(hObject, eventdata, handles)
	val = read_ui_str(hObject, 'n');
	if isempty(val) || (val < 0)
		warning('Invalid Sample Rate');
		update_ui_str(hObject, handles.Values.Export.SampleRate);
	else
		handles.Values.Export.SampleRate = val;
	end
	guidata(hObject, handles);
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function exportResample_ctrl_Callback(hObject, eventdata, handles)
	val = read_ui_val(hObject);
	if val
		enable_ui(handles.ExportSampleRate_label);
		enable_ui(handles.exportSampleRate_ctrl);
	else
		disable_ui(handles.ExportSampleRate_label);
		disable_ui(handles.exportSampleRate_ctrl);
	end
	handles.Values.Export.ResampleData = val;
	handles.Values.Export.SampleRate = ...
						read_ui_str(handles.exportSampleRate_ctrl, 'n');
	guidata(hObject, handles);
%-------------------------------------------------------------------------- 
%-------------------------------------------------------------------------- 
function exportFormat_ctrl_Callback(hObject, eventdata, handles)
	val = read_ui_val(hObject);
	switch(val)
		case 1
			handles.Values.Export.ExportFormat = 'MAT';
			fprintf('Export as MAT\n');
		case 2
			handles.Values.Export.ExportFormat = 'CSV';
			fprintf('Export as CSV\n');
	end
	guidata(hObject, handles);
%-------------------------------------------------------------------------- 
%-------------------------------------------------------------------------- 
%**************************************************************************
%**************************************************************************
%**************************************************************************


%**************************************************************************
%**************************************************************************
%**************************************************************************
% PLOT callbacks
%**************************************************************************
%**************************************************************************
%**************************************************************************
%-------------------------------------------------------------------------- 
%-------------------------------------------------------------------------- 
function Tmin_ctrl_Callback(hObject, eventdata, handles) %#ok<*DEFNU>
	newval = str2double(get(hObject, 'String'));
	if isnan(newval)
		 update_ui_str(hObject, handles.Values.Tmin);
		 errordlg('Tmin must be a number','Error');
	elseif newval > handles.Values.Tmax
		 update_ui_str(hObject, handles.Values.Tmin);
		 eStr = sprintf('Tmin must be less than Tmax (%.2f)', ...
									handles.Values.Tmax);
		 errordlg(eStr,'Error');
	else
		% Save the new Tmin_ctrl value
		handles.Values.Tmin = newval;
		guidata(hObject,handles)
		update_gui(hObject, handles);
	end
%-------------------------------------------------------------------------- 
%-------------------------------------------------------------------------- 
function Tmax_ctrl_Callback(hObject, eventdata, handles)
	newval = str2double(get(hObject, 'String'));
	if isnan(newval)
		 update_ui_str(hObject, handles.Values.Tmax);
		 errordlg('Tmax must be a number','Error');
	elseif newval < handles.Values.Tmin
		 update_ui_str(hObject, handles.Values.Tmax);
		 eStr = sprintf('Tmax must be greater than Tmin (%.2f)', ...
									handles.Values.Tmin);
		 errordlg(eStr,'Error');
	else
		% Save the new Tmax_ctrl value
		handles.Values.Tmax = newval;
		guidata(hObject,handles)
		update_gui(hObject, handles);
	end
%-------------------------------------------------------------------------- 
%-------------------------------------------------------------------------- 
function Ymin_ctrl_Callback(hObject, eventdata, handles)
	newval = str2double(get(hObject, 'String'));
	if isnan(newval)
		 update_ui_str(hObject, handles.Values.Ymin);
		 errordlg('Ymin must be a number','Error');
	elseif newval > handles.Values.Ymax
		 update_ui_str(hObject, handles.Values.Ymin);
		 eStr = sprintf('Ymin must be less than Ymax (%.2f)', ...
												handles.Values.Ymax);
		 errordlg(eStr,'Error');
	else
		% Save the new Ymin_ctrl value
		handles.Values.Ymin = newval;
		guidata(hObject,handles)
		update_gui(hObject, handles);
	end
%-------------------------------------------------------------------------- 
%-------------------------------------------------------------------------- 
function Ymax_ctrl_Callback(hObject, eventdata, handles)
	newval = str2double(get(hObject, 'String'));
	if isnan(newval)
		 update_ui_str(hObject, handles.Values.Ymax);
		 errordlg('Ymax must be a number','Error');
	elseif newval < handles.Values.Ymin
		 update_ui_str(hObject, handles.Values.Ymax);
		 eStr = sprintf('Ymax must be greater than Ymin (%.2f)', ...
									handles.Values.Ymin);
		 errordlg(eStr,'Error');
	else
		% Save the new Ymax_ctrl value
		handles.Values.Ymax = newval;
		guidata(hObject,handles)
		update_gui(hObject, handles);
	end
%-------------------------------------------------------------------------- 
%-------------------------------------------------------------------------- 
function Decimate_ctrl_Callback(hObject, eventdata, handles)
	% read value from ctrl string
	newval = read_ui_str(hObject, 'n');
	% check it!
	if newval < 1
		warning(['%s: Decimation Factor must be an ' ...
						'integer greater than 0'], mfilename);
		update_ui_str(hObject, handles.Values.Decifactor);
		return
	end
	handles.Values.Decifactor =  newval;
	guidata(hObject, handles);
	if handles.E.isInitialized
		update_gui(hObject, handles);
	end
%-------------------------------------------------------------------------- 
%-------------------------------------------------------------------------- 
function Grid_ctrl_Callback(hObject, eventdata, handles)
%----------------------------------------------------
% show or hide grid on plot
%----------------------------------------------------
	handles.Values.Grid = read_ui_val(hObject);
	guidata(hObject, handles);
	if handles.E.isInitialized
		update_gui(hObject, handles);
	end
%-------------------------------------------------------------------------- 
%-------------------------------------------------------------------------- 
function MousePos_ctrl_Callback(hObject, eventdata, handles)
%----------------------------------------------------
% turn on/off mouse position detect
%----------------------------------------------------
mpVal = read_ui_val(hObject);
	while mpVal == 1
		[ms, mv, button] = ginput(1);	
		update_ui_str(handles.Ms_text, sprintf('%.2f', ms));
		update_ui_str(handles.Mv_text, sprintf('%.2f', mv));
		if button == 3
			mpVal = 0;
			update_ui_val(hObject, 0);
		else
			mpVal = read_ui_val(hObject);
		end		
	end
%-------------------------------------------------------------------------- 
%-------------------------------------------------------------------------- 
function TraceMean_ctrl_Callback(hObject, eventdata, handles)
%----------------------------------------------------
% Executes when Show Mean checkbox is selected
%----------------------------------------------------
% Get value from object
	handles.Values.PlotMean = read_ui_val(hObject);
	sprintf('handles.Values.PlotMean = %f', handles.Values.PlotMean);
	guidata(hObject, handles);
	update_gui(hObject, handles);
%-------------------------------------------------------------------------- 
%-------------------------------------------------------------------------- 
function StimStartEnd_ctrl_Callback(hObject, eventdata, handles)
%----------------------------------------------------
% executes when Show Stim Onset/offset is selected
%----------------------------------------------------
	handles.Values.PlotStimStartEnd = read_ui_val(hObject);
	guidata(hObject, handles);
	if handles.E.isInitialized
		update_gui(hObject, handles);
	end
%-------------------------------------------------------------------------- 
%-------------------------------------------------------------------------- 
%**************************************************************************
%**************************************************************************
%**************************************************************************

%**************************************************************************
%**************************************************************************
%**************************************************************************
% --- Menu Callbacks
%**************************************************************************
%**************************************************************************
%**************************************************************************
%-------------------------------------------------------------------------- 
%-------------------------------------------------------------------------- 
function File_menu_Callback(hObject, eventdata, handles)
	LoadData(hObject, eventdata, handles);
%-------------------------------------------------------------------------- 
%-------------------------------------------------------------------------- 
function Plot_menu_Callback(hObject, eventdata, handles)
%-------------------------------------------------------------------------- 
%-------------------------------------------------------------------------- 
function Plot_Load_menu_Callback(hObject, eventdata, handles)
%-------------------------------------------------------------------------- 
%-------------------------------------------------------------------------- 
function Plot_Save_menu_Callback(hObject, eventdata, handles)
% 	[filename, pathname] = uiputfile('*.fig', 'Save plot as', handles.Files.rawpath);
	tmpfile = [handles.Files.basename '-C' ...
						num2str(handles.Values.CurrentCondition) '.fig'];
	[filename, pathname] = uiputfile('*.fig', 'Save plot as', ...
													fullfile(pwd, tmpfile));
	if isequal(filename, 0) || isequal(pathname, 0)
		disp('Cancelled Save Plot...')
		return
	end
	
	if handles.E.isInitialized
		h = figure;
		handles.E.PlotTracesForCondition(handles.Values.CurrentCondition, h);
		saveas(h, fullfile(pathname, filename), 'fig');
		close(h);
	end
%-------------------------------------------------------------------------- 
%*****************************************************************************
%*****************************************************************************
%*****************************************************************************


%**************************************************************************
%**************************************************************************
%**************************************************************************
% BUTTONS
%**************************************************************************
%**************************************************************************
%**************************************************************************
%-------------------------------------------------------------------------- 
% --- Executes on button press in calculate. (NOT used?)
%-------------------------------------------------------------------------- 
function calculate_Callback(hObject, eventdata, handles)
%-------------------------------------------------------------------------- 
%-------------------------------------------------------------------------- 
% --- Executes on button press in reset.
%-------------------------------------------------------------------------- 
function reset_Callback(hObject, eventdata, handles)
	initialize_gui(gcbf, handles, true);
%-------------------------------------------------------------------------- 
%-------------------------------------------------------------------------- 
% --- Executes on button press in debug.
function debug_Callback(hObject, eventdata, handles)
	keyboard
%-------------------------------------------------------------------------- 
%**************************************************************************
%**************************************************************************
%**************************************************************************

%*****************************************************************************
%*****************************************************************************
%*****************************************************************************
%*****************************************************************************
%*****************************************************************************
% --- Executes during object creation, after setting all properties.
%*****************************************************************************
%*****************************************************************************
%*****************************************************************************
function Tmin_ctrl_CreateFcn(hObject, eventdata, handles) %#ok<*INUSD>
	if ispc && isequal(get(hObject,'BackgroundColor'), ...
								get(0,'defaultUicontrolBackgroundColor'))
		 set(hObject,'BackgroundColor','white');
	end
function Tmax_ctrl_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), ...
								get(0,'defaultUicontrolBackgroundColor'))
		 set(hObject,'BackgroundColor','white');
	end
function Ymin_ctrl_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), ...
								get(0,'defaultUicontrolBackgroundColor'))
		 set(hObject,'BackgroundColor','white');
	end
function Ymax_ctrl_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), ...
								get(0,'defaultUicontrolBackgroundColor'))
		 set(hObject,'BackgroundColor','white');
	end
function Condition_ctrl_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), ...
								get(0,'defaultUicontrolBackgroundColor'))
		 set(hObject,'BackgroundColor','white');
	end
function Sweep_ctrl_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), ...
								get(0,'defaultUicontrolBackgroundColor'))
		 set(hObject,'BackgroundColor','white');
	end
function exportSampleRate_ctrl_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), ...
								get(0,'defaultUicontrolBackgroundColor'))
		 set(hObject,'BackgroundColor','white');
	end
function exportFormat_ctrl_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), ...
								get(0,'defaultUicontrolBackgroundColor'))
		 set(hObject,'BackgroundColor','white');
	end
function SpikeThreshold_ctrl_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), ...
								get(0,'defaultUicontrolBackgroundColor'))
		 set(hObject,'BackgroundColor','white');
	end
function SpikeHoldoff_ctrl_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), ...
								get(0,'defaultUicontrolBackgroundColor'))
		 set(hObject,'BackgroundColor','white');
	end
%-------------------------------------------------------------------------- 
%*****************************************************************************
%*****************************************************************************
