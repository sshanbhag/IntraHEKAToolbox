function [pul_tr, pul_rr] = readPUL(pulfile)
%-----------------------------------------------------------------------------
% [pul_tr, pul_rr] = readPUL(pulfile)
%-----------------------------------------------------------------------------
% 
% Reads data from HEKA .pul file
% 
%-----------------------------------------------------------------------------
% Input Arguments:
% 	pulfile		path/filename
% 
% Output Arguments:
% 	pul			pul data structure
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
MACRO_NAME_SIZE = 16; %#ok<NASGU>
COMMENT_SIZE = 80;
EPC9_STATS_SIZE = 104; %#ok<NASGU>
ROOT_TEXT_SIZE = 400;
BYTE = 'uint8';

%-----------------------------------------------------------------------------
%-----------------------------------------------------------------------------
% Read information from file
%-----------------------------------------------------------------------------
%-----------------------------------------------------------------------------
% open file
fp = fopen(pulfile, 'r', 'ieee-le');

% check if error opening file
if fp == -1
	error('%s: could not open file %s for reading', mfilename, pgffile);
end

% scan tree
pul_tr.magic_number = fread(fp, 1, 'int32');
pul_tr.number_tree_levels = fread(fp, 1, 'int32');

pul_tr.tree_level = zeros(pul_tr.number_tree_levels, 1);
for i = 1:pul_tr.number_tree_levels
  pul_tr.tree_level(i) = fread(fp, 1, 'int32');
end

% read root_record
pul_rr.version = fread(fp, 1, 'int16');
pul_rr.version_name14 = deblank(char(fread(fp, STR_SIZE, 'char')'));
pul_rr.file_name14 = deblank(char(fread(fp, STR_SIZE, 'char')'));
pul_rr.comments400 = deblank(char(fread(fp, ROOT_TEXT_SIZE, 'char')'));
pul_rr.start_time = fread(fp, 1, 'double');

pul_rr.number_children = fread(fp, 1, 'int32');   % number of children   -   groups

% read group records
for i = 1:pul_rr.number_children
	pul_rr.gr(i).label14 = deblank(char(fread(fp, STR_SIZE, 'char')'));
	pul_rr.gr(i).text80 = deblank(char(fread(fp, COMMENT_SIZE, 'char')'));
	pul_rr.gr(i).experimenter_number = fread(fp, 1, 'int32');
	pul_rr.gr(i).extra_long_real = fread(fp, 1, 'double');

	pul_rr.gr(i).number_children = fread(fp, 1, 'int32');     % number of children   -   series

	% read series records
	for j = 1:pul_rr.gr(i).number_children
		size_of = 0;
		pul_rr.gr(i).srr(j).time = fread(fp, 1, 'double');
		size_of = size_of + 8;
		pul_rr.gr(i).srr(j).bandwidth = fread(fp, 1, 'double');
		size_of = size_of + 8;
		pul_rr.gr(i).srr(j).pipette_potential = fread(fp, 1, 'double');
		size_of = size_of + 8;
		pul_rr.gr(i).srr(j).cell_potential = fread(fp, 1, 'double');
		size_of = size_of + 8;
		pul_rr.gr(i).srr(j).pipette_resistance = fread(fp, 1, 'double');
		size_of = size_of + 8;
		pul_rr.gr(i).srr(j).seal_resistance = fread(fp, 1, 'double');
		size_of = size_of + 8;
		pul_rr.gr(i).srr(j).background_noise = fread(fp, 1, 'double');
		size_of = size_of + 8;
		pul_rr.gr(i).srr(j).temperature = fread(fp, 1, 'double');
		size_of = size_of + 8;
		pul_rr.gr(i).srr(j).pipette_pressure = fread(fp, 1, 'double');
		size_of = size_of + 8;
		pul_rr.gr(i).srr(j).user_param1_value = fread(fp, 1, 'double');
		size_of = size_of + 8;
		pul_rr.gr(i).srr(j).user_param1_name14 = char(fread(fp, STR_SIZE, 'char')');
		size_of = size_of + STR_SIZE;
		pul_rr.gr(i).srr(j).user_param1_unit = fread(fp, 1, 'int16');
		size_of = size_of + 2;
		pul_rr.gr(i).srr(j).user_param2_value = fread(fp, 1, 'double');
		size_of = size_of + 8;
		pul_rr.gr(i).srr(j).user_param2_name14 = char(fread(fp, STR_SIZE, 'char')');
		size_of = size_of + STR_SIZE;
		pul_rr.gr(i).srr(j).user_param2_unit = fread(fp, 1, 'int16');
		size_of = size_of + 2;
		pul_rr.gr(i).srr(j).recording_mode = fread(fp, 1, BYTE);
		size_of = size_of + 1;
		pul_rr.gr(i).srr(j).filler1 = fread(fp, 1, BYTE);
		size_of = size_of + 1;
		pul_rr.gr(i).srr(j).comment80 = char(fread(fp, COMMENT_SIZE, 'char')');
		size_of = size_of + COMMENT_SIZE;
		% pul_rr.gr(i).srr(j).epc9_state104 = ReadString(b_reader, EPC9_STATS_SIZE)
		% size_of += EPC9_STATS_SIZE
		pul_rr.gr(i).srr(j).internal_solution = fread(fp, 1, 'int32');
		size_of = size_of + 4;
		pul_rr.gr(i).srr(j).external_solution = fread(fp, 1, 'int32');
		size_of = size_of + 4;
		pul_rr.gr(i).srr(j).extra_y_unit1 = fread(fp, 1, 'int16');
		size_of = size_of + 2;
		pul_rr.gr(i).srr(j).extra_y_unit2 = fread(fp, 1, 'int16');
		size_of = size_of + 2;
		pul_rr.gr(i).srr(j).disp_y_unit1 = fread(fp, 1, 'int32');
		size_of = size_of + 4;
		pul_rr.gr(i).srr(j).disp_y_unit2 = fread(fp, 1, 'int32');
		size_of = size_of + 4;
		pul_rr.gr(i).srr(j).fura_k = fread(fp, 1, 'double');
		size_of = size_of + 8;
		pul_rr.gr(i).srr(j).fura_min = fread(fp, 1, 'double');
		size_of = size_of + 8;
		pul_rr.gr(i).srr(j).fura_max = fread(fp, 1, 'double');
		size_of = size_of + 8;
		pul_rr.gr(i).srr(j).lock_in_ext_phase = fread(fp, 1, 'double');
		size_of = size_of + 8;
		pul_rr.gr(i).srr(j).timer = fread(fp, 1, 'double');
		size_of = size_of + 8;
		pul_rr.gr(i).srr(j).extra_long_real = fread(fp, 1, 'double');
		size_of = size_of + 8;
		buffer = fread(fp, pul_tr.tree_level(3) - size_of, BYTE); %#ok<NASGU>
		% number of children   -   sweeps
		pul_rr.gr(i).srr(j).number_children = fread(fp, 1, 'int32');
		
		% read sweep records
		for k = 1:pul_rr.gr(i).srr(j).number_children
			size_of = 0;
			pul_rr.gr(i).srr(j).swr(k).time = fread(fp, 1, 'double');
			size_of = size_of + 8;
			pul_rr.gr(i).srr(j).swr(k).stim_count = fread(fp, 1, 'int32');
			size_of = size_of + 4;
			pul_rr.gr(i).srr(j).swr(k).sweep_count = fread(fp, 1, 'int32');
			size_of = size_of + 4;
			pul_rr.gr(i).srr(j).swr(k).average_count = fread(fp, 1, 'int32');
			size_of = size_of + 4;
			pul_rr.gr(i).srr(j).swr(k).leak = fread(fp, 1, BYTE);
			size_of = size_of + 1;
			pul_rr.gr(i).srr(j).swr(k).second_trace = fread(fp, 1, BYTE);
			size_of = size_of + 1;
			pul_rr.gr(i).srr(j).swr(k).label14 = char(fread(fp, STR_SIZE, 'char')');
			size_of = size_of + STR_SIZE;
			pul_rr.gr(i).srr(j).swr(k).data_points = fread(fp, 1, 'int32');
			size_of = size_of + 4;
			pul_rr.gr(i).srr(j).swr(k).data = fread(fp, 1, 'int32');
			size_of = size_of + 4;
			pul_rr.gr(i).srr(j).swr(k).data_pointer = fread(fp, 1, 'int32');
			size_of = size_of + 4;
			pul_rr.gr(i).srr(j).swr(k).data_factor1 = fread(fp, 1, 'double');
			size_of = size_of + 8;
			pul_rr.gr(i).srr(j).swr(k).data_factor2 = fread(fp, 1, 'double');
			size_of = size_of + 8;
			pul_rr.gr(i).srr(j).swr(k).c_slow = fread(fp, 1, 'double');
			size_of = size_of + 8;
			pul_rr.gr(i).srr(j).swr(k).g_series = fread(fp, 1, 'double');
			size_of = size_of + 8;
			pul_rr.gr(i).srr(j).swr(k).rs_value = fread(fp, 1, 'double');
			size_of = size_of + 8;
			pul_rr.gr(i).srr(j).swr(k).rs_fraction = fread(fp, 1, 'double');
			size_of = size_of + 8;
			pul_rr.gr(i).srr(j).swr(k).zero_current = fread(fp, 1, 'double');
			size_of = size_of + 8;
			pul_rr.gr(i).srr(j).swr(k).online_y_result = fread(fp, 1, 'double');
			size_of = size_of + 8;
			pul_rr.gr(i).srr(j).swr(k).online_x_result = fread(fp, 1, 'double');
			size_of = size_of + 8;
			pul_rr.gr(i).srr(j).swr(k).total_points = fread(fp, 1, 'int32');
			size_of = size_of + 4;
			pul_rr.gr(i).srr(j).swr(k).offset = fread(fp, 1, 'int32');
			size_of = size_of + 4;
			pul_rr.gr(i).srr(j).swr(k).sweep_kind = fread(fp, 1, 'int16');
			size_of = size_of + 2;
			pul_rr.gr(i).srr(j).swr(k).fura_points = fread(fp, 1, 'int32');
			size_of = size_of + 4;
			pul_rr.gr(i).srr(j).swr(k).fura_data = fread(fp, 1, 'int32');
			size_of = size_of + 4;
			pul_rr.gr(i).srr(j).swr(k).fura_pointer = fread(fp, 1, 'int32');
			size_of = size_of + 4;
			pul_rr.gr(i).srr(j).swr(k).online_y_result2 = fread(fp, 1, 'double');
			size_of = size_of + 8;
			pul_rr.gr(i).srr(j).swr(k).online_x_result2 = fread(fp, 1, 'double');
			size_of = size_of + 8;
			pul_rr.gr(i).srr(j).swr(k).disp_factor1 = fread(fp, 1, 'double');
			size_of = size_of + 8;
			pul_rr.gr(i).srr(j).swr(k).disp_factor2 = fread(fp, 1, 'double');
			size_of = size_of + 8;
			pul_rr.gr(i).srr(j).swr(k).data_format = fread(fp, 1, BYTE);
			size_of = size_of + 1;
			pul_rr.gr(i).srr(j).swr(k).data_abscissa = fread(fp, 1, BYTE);
			size_of = size_of + 1;
			pul_rr.gr(i).srr(j).swr(k).timer = fread(fp, 1, 'double');
			size_of = size_of + 8;
			pul_rr.gr(i).srr(j).swr(k).sweep_spares = fread(fp, 1, 'int32');
			size_of = size_of + 4;
			buffer = fread(fp, pul_tr.tree_level(4) - size_of, BYTE); %#ok<NASGU>
			% number of children   -     = 0  now
			pul_rr.gr(i).srr(j).swr(k).number_children = fread(fp, 1, 'int32');
		end
	end
end


fclose(fp);
