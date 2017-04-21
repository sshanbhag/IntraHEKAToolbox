function pgf = readPGF(pgffile)
%-----------------------------------------------------------------------------
% pgf = readPGF(pgffile)
%-----------------------------------------------------------------------------
% IntraHEKAToolbox
%-------------------------------------------------------------------------
% 
% Reads data from HEKA .pgf file
% 
%-----------------------------------------------------------------------------
% Input Arguments:
% 	pgffile		path/filename
% 
% Output Arguments:
% 	pgf			pgf data structure
%
%-----------------------------------------------------------------------------
% See also: 
%-----------------------------------------------------------------------------

%-----------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neomed.edu
%-----------------------------------------------------------------------------
% Created: 21 March, 2012 (SJS)
% 				based on Visual Basic code clsOData.vb by Olga Galazyuk
%
% Revisions:
%	9 Jun 2016 (SJS): updated comments, added check for failed fopen()
%-----------------------------------------------------------------------------
% TO DO:
%-----------------------------------------------------------------------------



%-----------------------------------------------------------------------------
%-----------------------------------------------------------------------------
% some definitions
%-----------------------------------------------------------------------------
%-----------------------------------------------------------------------------
STR_SIZE = 14;
MACRO_NAME_SIZE = 16;
BYTE = 'uint8';

%-----------------------------------------------------------------------------
%-----------------------------------------------------------------------------
% Read information from file
%-----------------------------------------------------------------------------
%-----------------------------------------------------------------------------

% open file
fp = fopen(pgffile, 'r', 'ieee-le');

% check if error opening file
if fp == -1
	error('%s: could not open file %s for reading', mfilename, pgffile);
end

pgf.magic_number = fread(fp, 1, 'int32');
pgf.number_tree_levels = fread(fp, 1, 'int32');
pgf.tree_level = zeros(pgf.number_tree_levels, 1);
for i = 1:pgf.number_tree_levels
	pgf.tree_level(i) = fread(fp, 1, 'int32');
end

pgf.root_version = fread(fp, 1, 'int16');
% number of  children - stimulations
pgf.number_children = fread(fp, 1, 'int32');


% read stimulation records
for i = 1:pgf.number_children
	size_of = 0;
	tmp = fread(fp, STR_SIZE, 'char');
	smr(i).file_name14 = char(tmp'); %#ok<*AGROW>
	size_of = size_of + STR_SIZE;
	tmp = fread(fp, STR_SIZE, 'char');
	smr(i).entry_name14 = char(tmp');
	size_of = size_of + STR_SIZE;
	% sampling rate period
	smr(i).sample_interval = fread(fp, 1, 'double');
	size_of = size_of + 8;
	smr(i).filter_factor = fread(fp, 1, 'double');
	size_of = size_of + 8;
	smr(i).st_sweep_interval = fread(fp, 1, 'double');
	size_of = size_of + 8;
	smr(i).number_sweeps = fread(fp, 1, 'int32');
	size_of = size_of + 4;
	smr(i).number_repeats = fread(fp, 1, 'int32');
	size_of = size_of + 4;
	smr(i).repeat_wait = fread(fp, 1, 'double');
	size_of = size_of + 8;
	tmp = fread(fp, STR_SIZE, 'char');
	smr(i).linked_sequence14 = char(tmp');
	size_of = size_of + STR_SIZE;
	smr(i).linked_wait = fread(fp, 1, 'double');
	size_of = size_of + 8;
	smr(i).leak_count = fread(fp, 1, 'int32');
	size_of = size_of + 4;
	smr(i).leak_size = fread(fp, 1, 'double');
	size_of = size_of + 8;
	smr(i).leak_holding = fread(fp, 1, 'double');
	size_of = size_of + 8;
	smr(i).leak_alternate = fread(fp, 1, BYTE);
	size_of = size_of + 1;
	smr(i).alt_leak_averaging = fread(fp, 1, BYTE);
	size_of = size_of + 1;
	smr(i).leak_delay = fread(fp, 1, 'double');
	size_of = size_of + 8;
	smr(i).trigger_segment1 = fread(fp, 1, 'int16');
	size_of = size_of + 2;
	smr(i).trigger_time1 = fread(fp, 1, 'double');
	size_of = size_of + 8;
	smr(i).trigger_length1 = fread(fp, 1, 'double');
	size_of = size_of + 8;
	smr(i).trigger_amplitude1 = fread(fp, 1, 'double');
	size_of = size_of + 8;
	smr(i).trigger_dac1 = fread(fp, 1, 'int16');
	size_of = size_of + 2;
	smr(i).trigger_segment2 = fread(fp, 1, 'int16');
	size_of = size_of + 2;
	smr(i).trigger_time2 = fread(fp, 1, 'double');
	size_of = size_of + 8;
	smr(i).trigger_length2 = fread(fp, 1, 'double');
	size_of = size_of + 8;
	smr(i).trigger_amplitude2 = fread(fp, 1, 'double');
	size_of = size_of + 8;
	smr(i).trigger_dac2 = fread(fp, 1, 'int16');
	size_of = size_of + 2;
	smr(i).trigger_segment3 = fread(fp, 1, 'int16');
	size_of = size_of + 2;
	smr(i).trigger_time3 = fread(fp, 1, 'double');
	size_of = size_of + 8;
	smr(i).trigger_length3 = fread(fp, 1, 'double');
	size_of = size_of + 8;
	smr(i).trigger_amplitude3 = fread(fp, 1, 'double');
	size_of = size_of + 8;
	smr(i).trigger_dac3 = fread(fp, 1, 'int16');
	size_of = size_of + 2;
	smr(i).number_of_triggers = fread(fp, 1, 'int16');
	size_of = size_of + 2;
	smr(i).relevant_x_segment = fread(fp, 1, 'int16');
	size_of = size_of + 2;
	smr(i).relevant_y_segment = fread(fp, 1, 'int16');
	size_of = size_of + 2;
	smr(i).write_mode = fread(fp, 1, BYTE);
	size_of = size_of + 1;
	smr(i).increment_mode = fread(fp, 1, BYTE);
	size_of = size_of + 1;
	smr(i).total_sweep_length = fread(fp, 1, 'int32');
	size_of = size_of + 4;
	smr(i).max_sweep_length = fread(fp, 1, 'int32');
	size_of = size_of + 4;
	smr(i).input_channels = fread(fp, 1, 'int16');
	size_of = size_of + 2;
	smr(i).g_update = fread(fp, 1, BYTE);
	size_of = size_of + 1;
	smr(i).rel_abs_pot = fread(fp, 1, BYTE);
	size_of = size_of + 1;
	smr(i).has_continuous = fread(fp, 1, BYTE);
	size_of = size_of + 1;
	smr(i).log_increment = fread(fp, 1, BYTE);
	size_of = size_of + 1;
	smr(i).pgf_stim_dac = fread(fp, 1, 'int16');
	size_of = size_of + 2;
	smr(i).adc1 = fread(fp, 1, 'int16');
	size_of = size_of + 2;
	smr(i).adc2 = fread(fp, 1, 'int16');
	size_of = size_of + 2;
	smr(i).y_unit1 = fread(fp, 1, 'int16');
	size_of = size_of + 2;
	smr(i).y_unit2 = fread(fp, 1, 'int16');
	size_of = size_of + 2;
	smr(i).v_memb_increment = fread(fp, 1, 'int32');
	size_of = size_of + 4;
	smr(i).ext_trigger = fread(fp, 1, BYTE);
	size_of = size_of + 1;
	smr(i).file_template = fread(fp, 1, BYTE);
	size_of = size_of + 1;
	smr(i).stim_kind = fread(fp, 1, 'int16');
	size_of = size_of + 2;
	smr(i).lock_in_cycle = fread(fp, 1, 'double');
	size_of = size_of + 8;
	smr(i).lock_in_amplitude = fread(fp, 1, 'double');
	size_of = size_of + 8;
	smr(i).fura_on = fread(fp, 1, BYTE);
	size_of = size_of + 1;
	smr(i).v_memb_mode = fread(fp, 1, BYTE);     %pgf_st_vmemb_mode
	size_of = size_of + 1;
	smr(i).fura_tot_length = fread(fp, 1, 'double');
	size_of = size_of + 8;
	smr(i).fura_delay = fread(fp, 1, 'double');
	size_of = size_of + 8;
	smr(i).fura_length1 = fread(fp, 1, 'double');
	size_of = size_of + 8;
	smr(i).fura_length2 = fread(fp, 1, 'double');
	size_of = size_of + 8;
	smr(i).fura_wave_length0 = fread(fp, 1, 'double');
	size_of = size_of + 8;
	smr(i).fura_wave_length1 = fread(fp, 1, 'double');
	size_of = size_of + 8;
	smr(i).fura_wave_length2 = fread(fp, 1, 'double');
	size_of = size_of + 8;
	smr(i).fura_repeats = fread(fp, 1, 'int16');
	size_of = size_of + 2;
	smr(i).lock_in_skip = fread(fp, 1, 'int32');
	size_of = size_of + 4;
	smr(i).lock_in_v_reversal = fread(fp, 1, 'double');
	size_of = size_of + 8;
	smr(i).lock_in_mode = fread(fp, 1, BYTE);
	size_of = size_of + 1;
	smr(i).lock_in_show = fread(fp, 1, BYTE);
	size_of = size_of + 1;
	tmp = fread(fp, MACRO_NAME_SIZE, 'char');
	smr(i).config_macro16 = char(tmp');
	size_of = size_of + MACRO_NAME_SIZE;
	tmp = fread(fp, MACRO_NAME_SIZE, 'char');
	smr(i).end_macro16 = char(tmp');
	size_of = size_of + MACRO_NAME_SIZE;
	smr(i).ampl_mode_kind = fread(fp, 1, BYTE);
	size_of = size_of + 1;
	smr(i).filler1 = fread(fp, 1, BYTE);
	size_of = size_of + 1;
	
	buffer = fread(fp, pgf.tree_level(2) - size_of, BYTE); %#ok<NASGU> % ????????????????????

	smr(i).number_children = fread(fp, 1, 'int32');   % number of children - segments
	% read stimulation records
	for j = 1:smr(i).number_children
		size_of = 0;
		smr(i).sgr(j).se_class = fread(fp, 1, BYTE);
		size_of = size_of + 1;
		smr(i).sgr(j).is_holding = fread(fp, 1, BYTE);
		size_of = size_of + 1;
		smr(i).sgr(j).voltage = fread(fp, 1, 'double');
		size_of = size_of + 8;
		smr(i).sgr(j).duration = fread(fp, 1, 'double');
		size_of = size_of + 8;
		smr(i).sgr(j).delta_v_factor = fread(fp, 1, 'double');
		size_of = size_of + 8;
		smr(i).sgr(j).delta_v_increment = fread(fp, 1, 'double');
		size_of = size_of + 8;
		smr(i).sgr(j).delta_t_factor = fread(fp, 1, 'double');
		size_of = size_of + 8;
		smr(i).sgr(j).delta_t_increment = fread(fp, 1, 'double');
		size_of = size_of + 8;
		buffer = fread(fp, pgf.tree_level(3) - size_of, BYTE); %#ok<NASGU>
		smr(i).sgr(j).number_children = fread(fp, 1, 'int32');
	end

end


fclose(fp);

pgf.smr = smr;

