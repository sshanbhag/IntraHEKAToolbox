%-----------------------------------------------------------------------------
% statistic
%-----------------------------------------------------------------------------
% IntraHEKA Toolbox
% Class Definition
%-----------------------------------------------------------------------------
% Abstract statistic class
%-----------------------------------------------------------------------------
% See also: 
%-----------------------------------------------------------------------------

%-----------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neomed.edu
%-----------------------------------------------------------------------------
% Created: 4 May, 2012 (SJS)
%
% Revisions:
%-----------------------------------------------------------------------------
% TO DO:
%-----------------------------------------------------------------------------

classdef statistic
	properties
		Name
		Nvalues
		Values
	end
	
	methods
		%------------------------------------------------------------------------
		% Constructor
		%------------------------------------------------------------------------
		function s = statistic(varargin)
			nargs = size(varargin, 2);

			if nargs == 0
				s.Name = '';
				s.Nvalues = 0;
				s.Values = [];
				return
			end
			
			if nargs ~= 2
				error('statistic.%s: improper inputs', mfilename);
			end
			
			s.Name = varargin{1};
			s.Nvalues = size(varargin{2});
			s.Values = varargin{2};
		end
		%------------------------------------------------------------------------
		
		function SetValues(obj, NewValues)
			
			if isempty(NewValues)
				disp('empty')
				obj.Nvalues = 0;
				obj.Values = [];
				return
			end
						
			obj.Nvalues = size(NewValues);
			obj.Values = NewValues
		end
		
	end
end