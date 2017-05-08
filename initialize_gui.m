function initialize_gui(fig_handle, handles, isreset)
%------------------------------------------------------------------------
% initialize_gui(fig_handle, handles, isreset)
%------------------------------------------------------------------------
% IntraHEKAToolbox
%------------------------------------------------------------------------
% 
% function that initializes settings, values and GUI for HEKAview program
%
%------------------------------------------------------------------------
% Input Arguments:
%	fig_handle		handle to GUI figure
%	handles			handles struct from HEKAview
%	isreset			'0' for no reset, 1 for reset
% 
%------------------------------------------------------------------------
% See also: HEKAview
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: May 2017: separated from main HEKAview program file (SJS)
%
% Revisions:
%	8 Maym 2017 (SJS): added documentation, 
% 								code for OverlaySweeps, RemoveSpikesFromMean option
%------------------------------------------------------------------------
% TO DO:
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% If the Values field is present and the reset flag is false, it means
% we are we are just re-initializing a GUI by calling it from the cmd line
% while it is up. So, bail out as we dont want to reset the data.
%------------------------------------------------------------------------
if isfield(handles, 'Values') && ~isreset
	 return;
end

%------------------------------------------------------------------------
% set some default values
%------------------------------------------------------------------------
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

%------------------------------------------------------------------------
% experiment object
%------------------------------------------------------------------------
handles.E = experiment;

%------------------------------------------------------------------------
% line , text colors
%------------------------------------------------------------------------
handles.Settings.MeanLineColor = 'm';
handles.Settings.MeanTextColor = handles.Settings.MeanLineColor;
handles.Settings.StimStartLineColor = 'g';
handles.Settings.StimEndLineColor = 'r';

%------------------------------------------------------------------------
% set some initial values for plots
%------------------------------------------------------------------------
handles.Values.Tmin = 0;
handles.Values.Tmax = 1000;
handles.Values.Ymin = -100;
handles.Values.Ymax = 10;
handles.Values.Decifactor = 10;
handles.Values.Grid = 1;
handles.Values.PlotMean = 1;
handles.Values.PlotStimStartEnd = 1;
% these two values are used for plotting and computing mean
handles.Values.OverlaySweeps = 1;
update_ui_val(handles.OverlaySweeps_ctrl, handles.Values.OverlaySweeps);
handles.Values.RemoveSpikesFromMean = 1;
update_ui_val(handles.RemoveSpikesFromMean_ctrl, ...
					handles.Values.RemoveSpikesFromMean);

%------------------------------------------------------------------------
% spike detect settings
%------------------------------------------------------------------------
% threshold for spikes, mV
handles.Values.SpikeThreshold = -15;
update_ui_str(handles.SpikeThreshold_ctrl, handles.Values.SpikeThreshold);
% spike holdoff period, ms
handles.Values.HoldoffTime = 2;	
update_ui_str(handles.SpikeHoldoff_ctrl, handles.Values.HoldoffTime);

%------------------------------------------------------------------------
% export settings
%------------------------------------------------------------------------
handles.Values.Export.RemoveSpikes = 0;
handles.Values.Export.ResampleData = 0;
handles.Values.Export.SampleRate = [];
handles.Values.Export.ExportFormat = 'MAT';
update_ui_val(handles.exportFormat_ctrl, 1);
disable_ui(handles.exportSampleRate_ctrl);
disable_ui(handles.ExportSampleRate_label);

%------------------------------------------------------------------------
% store settings
%------------------------------------------------------------------------
guidata(handles.figure1, handles);

%------------------------------------------------------------------------
% update ctrls
%------------------------------------------------------------------------
update_ui_str(handles.Tmin_ctrl, handles.Values.Tmin);
update_ui_str(handles.Tmax_ctrl, handles.Values.Tmax);
update_ui_str(handles.Ymin_ctrl, handles.Values.Ymin);
update_ui_str(handles.Ymax_ctrl, handles.Values.Ymax);
update_ui_str(handles.Decimate_ctrl, handles.Values.Decifactor);
update_ui_val(handles.Grid_ctrl, handles.Values.Grid);
update_ui_val(handles.TraceMean_ctrl, handles.Values.PlotMean);
update_ui_val(handles.StimStartEnd_ctrl, ...
								handles.Values.PlotStimStartEnd);
update_ui_val(handles.exportRemoveSpikes_ctrl, ...
									handles.Values.Export.RemoveSpikes);
update_ui_val(handles.exportResample_ctrl, ...
									handles.Values.Export.ResampleData);
update_ui_str(handles.exportSampleRate_ctrl, ...
									handles.Values.Export.SampleRate);

%------------------------------------------------------------------------
% Update handles structure
%------------------------------------------------------------------------
guidata(handles.figure1, handles);
