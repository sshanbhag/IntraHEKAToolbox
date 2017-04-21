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

% Last Modified by GUIDE v2.5 17-Oct-2016 14:21:54

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

%-------------------------------------------------------------------------- 
%-------------------------------------------------------------------------- 
%-------------------------------------------------------------------------- 
function initialize_gui(fig_handle, handles, isreset)
	% If the metricdata field is present and the reset flag is false, it means
	% we are we are just re-initializing a GUI by calling it from the cmd line
	% while it is up. So, bail out as we dont want to reset the data.
	if isfield(handles, 'Values') && ~isreset
		 return;
	end
	% set some default values
	handles.Files.logpath = ...
			'/Users/sshanbhag/Work/Data/Mouse/IntracellularAmygdala/RawData';
	handles.Files.logfile = 'IntracellDataLog.csv';
	handles.Files.logname = fullfile(handles.Files.logpath, ...
																handles.Files.logfile);
	set(handles.textLogFile, 'String', handles.Files.logname);
	handles.Files.rawpath = ...
		'/Users/sshanbhag/Work/Data/IntracellularAmygdala/RawData/2012-03-12';
	handles.Files.basename = '2012-03-12-4'; 
	handles.Files.datname = [handles.Files.basename '.dat'];
	handles.Files.pgfname = [handles.Files.basename '.pgf'];
	handles.Files.pulname = [handles.Files.basename '.pul'];
	handles.Files.dscname = [handles.Files.basename '.dsc'];
	% experiment object
	handles.E = experiment;
	% line , text colors
	handles.Settings.MeanLineColor = 'm';
	handles.Settings.MeanTextColor = handles.Settings.MeanLineColor;
	handles.Settings.StimStartLineColor = 'g';
	handles.Settings.StimEndLineColor = 'r';
	% set some initial values
	handles.Values.Tmin = 0;
	handles.Values.Tmax = 1000;
	handles.Values.Ymin = -100;
	handles.Values.Ymax = 10;
	handles.Values.Decifactor = 10;
	handles.Values.Grid = 1;
	handles.Values.PlotMean = 1;
	handles.Values.PlotStimStartEnd = 1;
	% store settings
	guidata(handles.figure1, handles);
	% update ctrls
	update_ui_str(handles.Tmin_ctrl, handles.Values.Tmin);
	update_ui_str(handles.Tmax_ctrl, handles.Values.Tmax);
	update_ui_str(handles.Ymin_ctrl, handles.Values.Ymin);
	update_ui_str(handles.Ymax_ctrl, handles.Values.Ymax);
	update_ui_str(handles.Decimate_ctrl, handles.Values.Decifactor);
	update_ui_val(handles.Grid_ctrl, handles.Values.Grid);
	update_ui_val(handles.TraceMean_ctrl, handles.Values.PlotMean);
	update_ui_val(handles.StimStartEnd_ctrl, ...
									handles.Values.PlotStimStartEnd);
	% Update handles structure
	guidata(handles.figure1, handles);
%-------------------------------------------------------------------------- 
%-------------------------------------------------------------------------- 

%-------------------------------------------------------------------------- 
%-------------------------------------------------------------------------- 
%-------------------------------------------------------------------------- 
function update_gui(hObj, handles)
	%----------------------------------------------------------
	% select main plot axes
	%----------------------------------------------------------
	axes(handles.axesMain);
	%----------------------------------------------------------
	% check if data obj is initialized and ready
	%----------------------------------------------------------
	if ~handles.E.isInitialized
		warning('%s: Not initialized!', mfilename);
		return
	else
		% make local copy of handles to make things a little more compact
		% when calling parametersw
		H = handles;
	end
	%----------------------------------------------------------
	% update information amount stimulus/condition/experiment
	%----------------------------------------------------------
	S = H.E.GetStimulusForCondition(H.Values.CurrentCondition);
	update_ui_str(H.Animal_text, H.E.Info.Animal);
	update_ui_str(H.Depth_text, H.E.Info.Depth);
	update_ui_str(H.AuditoryStim_text, S.AuditoryStimulus);
	update_ui_str(H.OtherStim_text, S.OtherStimulus);
	update_ui_str(H.Comments_text, S.Comments);	
	%----------------------------------------------------------
	% update plot depending on what is to be plotted...
	%----------------------------------------------------------
	if H.Values.CurrentSweep == -1
		% SHOW ALL TRACES
		% plot all traces for this condition
		H.E.PlotTracesForCondition(	H.Values.CurrentCondition, ...
														H.Values.Decifactor, ...
														gcf);
		% set x axis (time) limits
		xlim([H.Values.Tmin H.Values.Tmax]);
		% turn off grid
		grid off
	else
		% SHOW SINGLE TRACE
		% get trace data and plot them
		axes(H.axesMain);
		[dt, S, T] = H.E.GetResampledTrace(	H.Values.CurrentSweep, ...
																H.Values.Decifactor);
		tvec = 1000 .* ((1:length(S)) - 1)  .* dt;
		plot(tvec, H.Values.Ymax - (1.5 + normalize(S)), ...
															'k', tvec, 1000 * T, 'b');
		% set plot limits
		ylim([H.Values.Ymin H.Values.Ymax]);
		xlim([H.Values.Tmin H.Values.Tmax]);
		% enable grid if Values.Grid is set
		if H.Values.Grid
			grid on
			drawnow
			update_ui_val(H.Grid_ctrl, 1);
		else
			grid off
			drawnow
			update_ui_val(H.Grid_ctrl, 0);
		end
		% show line for mean of trace
		if H.Values.PlotMean
			[t1, t2] = H.E.GetTrace(H.Values.CurrentSweep); %#ok<ASGLU>
			% get avg and std of trace in mV (hence the factor of 1000)
			t2avg = 1000 * mean(t2);
			t2std = std(1000 * t2); %#ok<NASGU>
			%		axes(H.axesMain);
			% draw line
			line(xlim, [t2avg t2avg],  'Color', H.Settings.MeanLineColor);
			text(H.Values.Tmax, t2avg, sprintf(' %.2f', t2avg), 'Color', ...
															H.Settings.MeanTextColor);
		end
	end
	%----------------------------------------------------------
	% show stim start/end lines
	%----------------------------------------------------------
	if H.Values.PlotStimStartEnd
		stiminfo = H.E.GetStimParamForCondition(H.Values.CurrentCondition);
		% draw start line
		line(stiminfo.start .* [1 1], ylim, 'Color', ...
								H.Settings.StimStartLineColor);
		% draw end line
		line(stiminfo.end .* [1 1], ylim, 'Color', ...
								H.Settings.StimEndLineColor);
	end
%-------------------------------------------------------------------------- 
%-------------------------------------------------------------------------- 


%**************************************************************************
%**************************************************************************
%**************************************************************************
% CTRL callbacks
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

%-------------------------------------------------------------------------- 
%-------------------------------------------------------------------------- 
% --- Executes on Grid ctrl checkbox
%-------------------------------------------------------------------------- 
function Grid_ctrl_Callback(hObject, eventdata, handles)
	handles.Values.Grid = read_ui_val(hObject);
	guidata(hObject, handles);
	if handles.E.isInitialized
		update_gui(hObject, handles);
	end
%-------------------------------------------------------------------------- 

%-------------------------------------------------------------------------- 
%-------------------------------------------------------------------------- 
% --- Executes on Mouse Position button press
%-------------------------------------------------------------------------- 
function MousePos_ctrl_Callback(hObject, eventdata, handles)
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
%-------------------------------------------------------------------------- 
% Executes when Show Mean checkbox is selected
%-------------------------------------------------------------------------- 
function TraceMean_ctrl_Callback(hObject, eventdata, handles)
	% Get value from object
	handles.Values.PlotMean = read_ui_val(hObject);
	sprintf('handles.Values.PlotMean = %f', handles.Values.PlotMean);
	guidata(hObject, handles);
	update_gui(hObject, handles);
%-------------------------------------------------------------------------- 
%-------------------------------------------------------------------------- 

%-------------------------------------------------------------------------- 
%-------------------------------------------------------------------------- 
% executes when Show Stim Onset/offset is selected
%-------------------------------------------------------------------------- 
function StimStartEnd_ctrl_Callback(hObject, eventdata, handles)
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
% BUTTONS
%**************************************************************************
%**************************************************************************
%**************************************************************************

%-------------------------------------------------------------------------- 
% --- Executes on button press in calculate.
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

%-------------------------------------------------------------------------- 
% --- Executes on button press in buttonLoadData.
function buttonLoadData_Callback(hObject, eventdata, handles)
	LoadData(hObject, eventdata, handles);
%-------------------------------------------------------------------------- 

%-------------------------------------------------------------------------- 
% --- Executes on button press in ExportSweeps.
function ExportSweeps_Callback(hObject, eventdata, handles)
	%-----------------------------------------------
	% check init
	%-----------------------------------------------
	if ~handles.E.isInitialized
		warning('%s: Not initialized!', mfilename);
		return
	end
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
		disp('Cancelled Load File..')
		return
	end
	%-----------------------------------------------
	% get sweeps, info
	%-----------------------------------------------
	[Stimuli, Sweeps] = handles.E.GetTrace(handles.Values.Sweeplist); %#ok<ASGLU>
	Fs = handles.E.GetSampleRate; %#ok<NASGU>
	%-----------------------------------------------
	% export to MAT
	%-----------------------------------------------
	save(fname, 'Stimuli', 'Sweeps', 'Fs', '-MAT');
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% --- Executes on button press in buttonMean.
function buttonMean_Callback(hObject, eventdata, handles)
	%-----------------------------------------------
	% despike?
	%-----------------------------------------------
	yn = uiyesno('title', 'Sweep Mean', 'string', 'Remove Spikes?');
	% get mean trace (without or with spikes)
	if strcmpi(yn, 'yes')
		disp('Removing Spikes...');
		[sweep_mean, sweep_std, stim] = ...
										handles.E.MeanTraceForCondition( ...
										handles.Values.CurrentCondition, ...
										'NoSpikes');
	else
		[sweep_mean, sweep_std, stim] = ...
										handles.E.MeanTraceForCondition( ...
										handles.Values.CurrentCondition);
	end
	%-----------------------------------------------
	% plot mean and error bars
	%-----------------------------------------------
	% get sampling rate
	Fs = handles.E.GetSampleRate;
	% time vector
	tvec = 1000 * ((1:length(sweep_mean)) - 1) * (1/Fs);
	% create new fig
	figure
	% plot mean
	plot(tvec, 1000*sweep_mean);
	xlabel('Time (ms)')
	ylabel('mV');
	fname = sprintf('%s_%d-%d.mat', ...
								handles.Files.basename, ...
								handles.Values.Sweeplist(1), ...
								handles.Values.Sweeplist(end) );
	title(fname, 'Interpreter', 'none');
	% plot error bars
	hold on
		plot(tvec, 1000*(sweep_mean + sweep_std), 'r');
		plot(tvec, 1000*(sweep_mean - sweep_std), 'r');
	hold off
	% plot stimulus
	yminmax = ylim;
	stim2 = yminmax(2)+0.1*(normalize(stim) - 1);
	hold on
		plot(tvec, stim2, 'Color', 0.5 * [1 1 1])
	hold off
	ylim([yminmax(1) max(stim2)]);
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


%-------------------------------------------------------------------------- 
%-------------------------------------------------------------------------- 
% --- Load Data
%-------------------------------------------------------------------------- 
function LoadData(hObject, eventdata, handles)
	%-----------------------------------------------
	% get file from user
	%-----------------------------------------------
	[filename, pathname] = uigetfile('*.dat', 'Pick a HEKA data file', ...
													handles.Files.rawpath);
	% if user cancelled, abort
	if isequal(filename, 0) || isequal(pathname, 0)
		disp('Cancelled Load File..')
		return
	end
	%-----------------------------------------------
	% Update files information
	%-----------------------------------------------
	[rawpath, basename, ~] = fileparts(fullfile(pathname, filename));
	handles.Files.rawpath = rawpath;
	handles.Files.basename = basename;
	handles.Files.datname = [handles.Files.basename '.dat'];
	handles.Files.pgfname = [handles.Files.basename '.pgf'];
	handles.Files.pulname = [handles.Files.basename '.pul'];
	handles.Files.dscname = [handles.Files.basename '.dsc'];
	%-----------------------------------------------
	% read in information from files
	%-----------------------------------------------
	% read info from pgf file
	pgf = readPGF(fullfile(handles.Files.rawpath, handles.Files.pgfname));
	% read info from pul file
	pul = struct('tr', [], 'rr', []);
	[pul.tr, pul.rr] = readPUL(fullfile(handles.Files.rawpath, ...
										handles.Files.pulname));
	% build dsc struct
	dsc = buildDSC(fullfile(handles.Files.rawpath, ...
									handles.Files.datname), pgf, pul.rr); %#ok<NASGU>
	% create experiment object
	handles.E = experiment(handles.Files.rawpath, ...
										handles.Files.datname, 'initialize');
	% load log file
	try
		logData = readLog(handles.Files.logname);
	catch errMsg
		sprintf('%s', errMsg.identifier)
		error(['%s: error reading log file, ' ...
					'check path or filename and re-load'], mfilename);
	end
	% set information in experiment object from logData
	handles.E.SetInfoFromLog(logData);
	% store changes
	guidata(hObject, handles);
	%-----------------------------------------------
	% update Condition list box
	%-----------------------------------------------
	handles.Values.Nconditions = handles.E.Info.Nconditions;
	handles.Values.CurrentCondition = 1;
	update_ui_val(handles.Condition_ctrl, handles.Values.CurrentCondition);
	cstr = cell(handles.Values.Nconditions, 1);
	for n = 1:handles.Values.Nconditions 
		cstr{n} = sprintf('%d', n);
	end
	set(handles.Condition_ctrl, 'String', cstr);
	%-----------------------------------------------
	% update Sweep list box
	%-----------------------------------------------
	% first get the list of sweeps for current condition (and # of sweeps)
	sweeplist = handles.E.GetSweepListForCondition(...
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
	%-----------------------------------------------
	% update File string in Recording Info panel
	%-----------------------------------------------
	update_ui_str(handles.File_text, handles.Files.datname);
	%-----------------------------------------------
	% Store changes, update gui
	%-----------------------------------------------	
	guidata(hObject, handles);
	update_gui(hObject, handles);
%-------------------------------------------------------------------------- 
%-------------------------------------------------------------------------- 



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
	if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
		 set(hObject,'BackgroundColor','white');
	end
	
function Tmax_ctrl_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
		 set(hObject,'BackgroundColor','white');
	end	
function Ymin_ctrl_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
		 set(hObject,'BackgroundColor','white');
	end
function Ymax_ctrl_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
		 set(hObject,'BackgroundColor','white');
	end
function Condition_ctrl_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
		 set(hObject,'BackgroundColor','white');
	end
function Sweep_ctrl_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
		 set(hObject,'BackgroundColor','white');
	end
%-------------------------------------------------------------------------- 
%*****************************************************************************
%*****************************************************************************
