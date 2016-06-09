function h = PlotIndividualTracesForCondition(obj, Condition, fHandle)
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% plot all traces for condition
%------------------------------------------------------------------------
	% make sure object is initialized
	if ~obj.isInitialized
		warning('%s: object not initialized', mfilename);
		if nargout
			h = [];
		end
		return

	elseif ~isfield(obj.Info, 'Nconditions')
		warning('%s: Nconditions not set', mfilename);
		fprintf('\t\tuse SetInfoFromLog method to set information from log file\n\n');
		if nargout
			h = [];
		end
		return

		% check that Condition is inbounds
	elseif ~between(Condition, 1, obj.Info.Nconditions)
		warning('%s: Condition %d not within bounds (Nconditions: %d)', ...
						mfilename, Condition, obj.Info.Nconditions);
		if nargout
			h = [];
		end
		return
	end

	if exist('fHandle', 'var')
		figure(fHandle);
	end

	% get the sweeps for the indicated Condition and plot them
	sinfo = obj.Info.Stimulus(Condition);		
	sweeplist = obj.GetSweepListForCondition(Condition);
	[t1, t2] = obj.GetTrace(sweeplist);
	htmp = plotIndividualTraces(	t2, 'Stimulus', ...
											t1{1}, 'SampleInterval', ...
											obj.Sweeps(sweeplist(1)).Rate	);
	% create title, and display it
	titlestr =	{	obj.BaseName, ...
						obj.Info.Animal, ...
						[sinfo.AuditoryStimulus ' ' sinfo.OtherStimulus] ...
					};
	title(titlestr)

	% return handle if it was requested
	if nargout
		h = htmp;
	end
% end
%------------------------------------------------------------------------
%------------------------------------------------------------------------
