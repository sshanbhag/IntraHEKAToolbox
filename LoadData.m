
function LoadData(hObject, eventdata, handles)
%-------------------------------------------------------------------------- 
%-------------------------------------------------------------------------- 
% --- Load Data
%-------------------------------------------------------------------------- 
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
if isempty(read_ui_str(handles.exportSampleRate_ctrl))
	update_ui_str(handles.exportSampleRate_ctrl, handles.E.GetSampleRate);
end
update_ui_str(handles.SampleRate_text, ...
							sprintf('%.0f', handles.E.GetSampleRate));
%-----------------------------------------------
% Store changes, update gui
%-----------------------------------------------	
guidata(hObject, handles);
update_gui(hObject, handles);
