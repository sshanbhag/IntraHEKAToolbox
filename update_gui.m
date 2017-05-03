
function update_gui(hObj, handles)
%-------------------------------------------------------------------------- 
%-------------------------------------------------------------------------- 

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
update_ui_str(handles.SampleRate_text, ...
							sprintf('%d', H.E.GetSampleRate));
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
	% enable grid if Values.Grid is set
	if H.Values.Grid
		grid(H.axesMain, 'on');
		drawnow
		update_ui_val(H.Grid_ctrl, 1);
	else
		grid(H.axesMain, 'off');
		drawnow
		update_ui_val(H.Grid_ctrl, 0);
	end
else
	% SHOW SINGLE TRACE
	% get trace data and plot them
	axes(H.axesMain);
	[dt, S, T] = H.E.GetDecimatedTrace(	H.Values.CurrentSweep, ...
															H.Values.Decifactor);
	tvec = 1000 .* ((1:length(S)) - 1)  .* dt;
	plot(tvec, H.Values.Ymax - (1.5 + normalize(S)), ...
														'k', tvec, 1000 * T, 'b');
	% set plot limits
	ylim([H.Values.Ymin H.Values.Ymax]);
	xlim([H.Values.Tmin H.Values.Tmax]);
	% enable grid if Values.Grid is set
	if H.Values.Grid
		grid(H.axesMain, 'on');
		drawnow
		update_ui_val(H.Grid_ctrl, 1);
	else
		grid(H.axesMain, 'off');
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
