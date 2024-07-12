%DIPGETPREF   Gets a DIPimage preference
%   V = DIPGETPREF('name') retrieves the value of the named DIPimage
%   preference.
%   V = DIPGETPREF retrieves the value of all preferences in a struct V.
%   V = DIPGETPREF('factory') retrieves the default values of all
%   preferences.
%
%   When leaving out the output variable, DIPGETPREF prints the values
%   of the preferences to the command window. In this case, boolean
%   values are represented as 'on' or 'off'.
%
%   The property names and values are described in the user manual.
%
%   See also: DIPSETPREF

% (c)2017, Cris Luengo.
% (c)1999-2014, Delft University of Technology.
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
%
%    http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

function value = dipgetpref(name)

if nargin<1
   name = 'list';
end
if ~ischar(name)
   error('Input argument should be a string.')
end
if strcmpi(name,'factory')
   name = 'defaults';
end
switch lower(name)
   case {'list','defaults'}
      prefs = dippreferences(name);
      if nargout == 0
         printvalues(prefs);
      else
         value = prefs;
      end
   otherwise
      pref = dippreferences('get',name);
      if nargout == 0
         printvalue(name,pref);
      else
         value = pref;
      end
end

function printvalues(prefs)
snames = sort(fieldnames(prefs));
for ii=1:length(snames)
   value = subsref(prefs,substruct('.',snames{ii}));
   if islogical(value)
      value = printboolean(value);
   elseif ~ischar(value)
      value = mat2str(value);
   end
   disp(['    ',snames{ii},' = ',value])
end
disp(' ')

function printvalue(name,value)
if islogical(value)
   value = printboolean(value);
elseif ~ischar(value)
   value = mat2str(value);
end
disp(['    ',name,' = ',value])
disp(' ')

function str = printboolean(val)
if val
   str = 'on';
else
   str = 'off';
end
