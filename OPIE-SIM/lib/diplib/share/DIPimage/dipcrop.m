%DIPCROP   Interactive image cropping
%   [B,C] = DIPCROP(H) allows the user to select a rectangular region in
%   the image show in figure with handle H. H defaults to the current figure.
%
%   The cropped image is returned in B. If it is a 3D image, it is not
%   possible to select in the 3rd dimension, so that all slices are
%   returned.
%
%   In C the selected box is returned in the form: [upper left corner; size]:
%   B = image(C(1,1) + (0:C(2,1)), C(1,2) + (0:C(2,2)))
%
%   It is still possible to use all the menus in the victim figure
%   window, but you won't be able to access any of the tools (like
%   zooming and testing). You can, however, step through a 3D
%   image by using the keyboard.
%
%   Note: If you feel the need to interrupt this function with Ctrl-C,
%   you will need to refresh the display (by re-displaying the image
%   or changing the 'Actions' state).
%
%   See also DIPSHOW, DIPGETIMAGE, DIPGETCOORDS, DIPPROFILE.

% (c)2018, Cris Luengo.
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

function [out,roi] = dipcrop(fig)

% Parse input
if nargin == 0
   fig = get(0,'CurrentFigure');
   if isempty(fig)
      error('No figure window open to do operation on')
   end
else % nargin == 1
   try
      fig = getfigh(fig);
   catch
      error('Argument must be a valid figure handle')
   end
end

tag = get(fig,'Tag');
if ~strncmp(tag,'DIP_Image',9)
   error('DIPCROP only works on images displayed using DIPSHOW')
end
ax = findobj(fig,'Type','axes');
if length(ax)~=1
   error('DIPCROP only works on images displayed using DIPSHOW')
end

% Store old settings
wbdF = get(fig,'WindowButtonDownFcn');
wbuP = get(fig,'WindowButtonUpFcn');
wbmF = get(fig,'WindowButtonMotionFcn');
bdF = get(fig,'ButtonDownFcn');
pscd = get(fig,'PointerShapeCData');
pshs = get(fig,'PointerShapeHotSpot');
ptr = get(fig,'pointer');

% Set new settings
figure(fig);
set(fig,'WindowButtonDownFcn','',...
        'WindowButtonUpFcn','',...
        'WindowButtonMotionFcn','',...
        'ButtonDownFcn','');
dipfig_setpointer(fig,'cross');

% Do your stuff
done = 0;
while ~done
   set(fig,'WindowButtonDownFcn',@dipcropButtonDownFcn);
   waitfor(fig,'WindowButtonDownFcn'); % The ButtonUp callback changes the callback.
                                       % This way, we also detect a change in state!
   if ~ishandle(fig)
      error('You closed the window! That wasn''t the deal!')
   end
   if ~strcmp(get(fig,'WindowButtonDownFcn'),'Done!')
      % The user just changed the state. Store the new settings and revert to our own...
      wbdF = get(fig,'WindowButtonDownFcn');
      wbuP = get(fig,'WindowButtonUpFcn');
      wbmF = get(fig,'WindowButtonMotionFcn');
      bdF = get(fig,'ButtonDownFcn');
      pscd = get(fig,'PointerShapeCData');
      pshs = get(fig,'PointerShapeHotSpot');
      ptr = get(fig,'pointer');
      set(fig,'WindowButtonDownFcn','',...
              'WindowButtonUpFcn','',...
              'WindowButtonMotionFcn','',...
              'ButtonDownFcn','');
      dipfig_setpointer(fig,'cross');
   else
      done = 1;
   end
end
udata = get(fig,'UserData');
if isfield(udata,'rectangle')
   rect = udata.rectangle;
   udata = rmfield(udata,'rectangle');
   set(fig,'UserData',[]);
   set(fig,'UserData',udata);
else
   error('Rectangle selection failed.')
end
pt = rect([1,2]);
sz = max([1,1],rect([3,4])) - 1;
udata = get(fig,'UserData');
if length(udata.imsize)==1
   out = out(pt(1)+(0:sz(1)));
else
   out = imagedisplay(udata.handle,'input');
   slicing = imagedisplay(udata.handle,'slicing');
   subs = cell(ndims(out),1);
   subs(:) = {':'};
   subs{slicing(1)} = pt(1)+(0:sz(1));
   subs{slicing(2)} = pt(2)+(0:sz(2));
   out = out(subs{:});
end
roi = [pt;sz];

% Restore old settings
set(fig,'WindowButtonDownFcn',wbdF,...
        'WindowButtonUpFcn',wbuP,...
        'WindowButtonMotionFcn',wbmF,...
        'ButtonDownFcn',bdF,...
        'PointerShapeCData',pscd,...
        'PointerShapeHotSpot',pshs,...
        'pointer',ptr);


%
% Callback function for mouse down in image
%
function dipcropButtonDownFcn(fig,~)
if strncmp(get(fig,'Tag'),'DIP_Image',9)
   udata = get(fig,'UserData');
   ax = findobj(fig,'Type','axes');
   if length(ax)~=1
      return
   end
   udata.ax = ax;
   udata.oldAxesUnits = get(ax,'Units');
   udata.oldNumberTitle = get(fig,'NumberTitle');
   set(ax,'Units','pixels');
   udata.coords = dipfig_getcurpos(ax); % Always over image!
   if length(udata.imsize)==1
      ylim = get(ax,'YLim');
      pos = [udata.coords(1)-0.5,ylim(1)-1,1,diff(ylim)+3];
   else
      pos = [udata.coords-0.5,1,1];
   end
   if useshg2
      udata.recth = rectangle('Position',pos,'EdgeColor',[0,0,0.8]);
   else
      udata.recth = rectangle('Position',pos,'EraseMode','xor','EdgeColor',[0,0,0.8]);
   end
   set(fig,'WindowButtonMotionFcn',@dipcropMotionFcn,...
           'WindowButtonUpFcn',@dipcropButtonUpFcn,...
           'NumberTitle','off',...
           'UserData',[]);   % Solve MATLAB bug!
   set(fig,'UserData',udata);
   updateDisplay(fig,ax,udata);
end


%
% Callback function for mouse move in image
%
function dipcropMotionFcn(fig,~)
if strncmp(get(fig,'Tag'),'DIP_Image',9)
   udata = get(fig,'UserData');
   updateDisplay(fig,udata.ax,udata);
end


%
% Callback function for mouse up in image
%
function dipcropButtonUpFcn(fig,~)
if strncmp(get(fig,'Tag'),'DIP_Image',9)
   udata = get(fig,'UserData');
   delete(udata.recth);
   pt = dipfig_getcurpos(udata.ax);
   udata.rectangle = [min(pt,udata.coords),abs(pt-udata.coords)+1];
   udata = rmfield(udata,{'recth','coords'});
   set(udata.ax,'Units',udata.oldAxesUnits);
   set(fig,'WindowButtonMotionFcn','','WindowButtonUpFcn','',...
           'NumberTitle',udata.oldNumberTitle);
   udata = rmfield(udata,{'ax','oldAxesUnits','oldNumberTitle'});
   dipfig_titlebar(fig,udata);
   set(fig,'UserData',[]);    % Solve MATLAB bug!
   set(fig,'UserData',udata,'WindowButtonDownFcn','Done!');
end


%
% Retrieve coordinates and image value, and update figure title.
%
function updateDisplay(fig,ax,udata)
pt = dipfig_getcurpos(ax);
if length(udata.imsize) == 1
   delta = abs(pt(1)-udata.coords(1))+1;
   pos = get(udata.recth,'Position');
   pos(1) = min(pt(1),udata.coords(1))-0.5;
   pos(3) = delta;
   set(udata.recth,'Position',pos);
   set(fig,'Name',['(',num2str(pt(1)),') ',' size: ',num2str(delta(1))]);
else
   delta = abs(pt-udata.coords)+1;
   set(udata.recth,'Position',[min(pt,udata.coords)-0.5,delta]);
   set(fig,'Name',['(',num2str(pt(1)),',',num2str(pt(2)),') ',...
          ' size: ',num2str(delta(1)),'x',num2str(delta(2))]);
end
