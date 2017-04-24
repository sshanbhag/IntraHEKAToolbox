%-----------------------------------------------------------------------------
% sweep
%-----------------------------------------------------------------------------
% IntraHEKA Toolbox
% Class Definition
%-----------------------------------------------------------------------------
% Some logic:
% 	each instance of sweep will know where its trace data are located within
% 	the binary .dat file, but will not store the trace data itself.  This is 
% 	to minimize the memory footprint of the data structure.
%-----------------------------------------------------------------------------
% See also: readPGF, readPUL, readHEKA, readDAT
%-----------------------------------------------------------------------------

%-----------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neomed.edu
%-----------------------------------------------------------------------------
% Created: 4 April, 2012 (SJS)
%
% Revisions:
%	12 Apr 2012 (SJS)
%		-	converted to value class (instead of handle - was causing 
%			unwanted behavior for arrays of sweep classes)
%-----------------------------------------------------------------------------
% TO DO:
%-----------------------------------------------------------------------------

classdef sweep
	properties
		Nsamples = [];
		Sweepindex = [];
		Rate = [];
		fp1 = [];
		fp2 = [];
		ID = [];
		IDname = '';
		Spiketimes = [];
	end
	
	methods
		%------------------------------------------------------------------------
		% Constructor
		%------------------------------------------------------------------------
		function s = sweep(nsamples, index, rate, fp1, fp2, varargin)
			optargin = size(varargin, 2);
			stdargin = nargin - optargin;

			if stdargin == 0
				s.Nsamples = [];
				s.Sweepindex = [];
				s.Rate = [];
				s.fp1 = [];
				s.fp2 = [];
				s.ID = [];
				s.IDname = '';
				s.Spiketimes = [];
				return
			end
				
			s.Nsamples = nsamples;
			s.Sweepindex = index;
			s.Rate = rate;
			s.fp1 = fp1;
			s.fp2 = fp2;
			if length(varargin) >= 1
				s.ID = varargin{1};
			end
			if length(varargin) > 1
				s.IDname = varargin{2};
			end
		end
		%------------------------------------------------------------------------
		
		%------------------------------------------------------------------------
		% read data
		%------------------------------------------------------------------------
		function [trace1, trace2] = ReadData(obj, datfile)
			trace1 = [];
			trace2 = [];

			if ~exist(datfile, 'file')
				warning('%s: data file %s does not exist', mfilename, datfile);
				return
			elseif any(isempty([obj.Nsamples obj.fp1 obj.fp2]))
				warning('%s: invalid or empty sweep information');
				return
			end
			
			try
				% open file
				fp = fopen(datfile, 'r', 'ieee-le');
				
				fseek(fp, obj.fp1, 'bof');
				trace1 = fread(fp, obj.Nsamples, 'int16');

				% read the data and scale using data_factor2 (??)
				fseek(fp, obj.fp2, 'bof');
				trace2 = fread(fp, obj.Nsamples, 'int16');
				
				fclose(fp);
			catch
				warning('%s: error while reading file %s', mfilename, datfile);
				trace1 = [];
				trace2 = [];
				return
			end
		end
		%------------------------------------------------------------------------
		
	end
end