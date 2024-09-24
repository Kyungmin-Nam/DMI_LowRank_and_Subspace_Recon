function varargout = clim(varargin)
%CLIM Set or query color limits
%   CLIM(limits) specifies the color limits for the current axes. Specify
%   limits as a two-element vector of the form [cmin cmax], where cmax is a
%   numeric value greater than cmin.
%
%   cl = CLIM returns a two-element vector containing the color limits for
%   the current axes.
%   
%   CLIM('auto') lets the axes choose the color limits. This command sets
%   the CLimMode property for the axes to 'auto'.
%
%   CLIM('manual') freezes the color limits at the current values.This 
%   command sets the XLimMode property for the axes to 'manual'.
% 
%   m = CLIM('mode') returns the current value of the color limits mode,
%   which is either 'auto' or 'manual'. By default, the mode is automatic
%   unless you specify limits or set the mode to manual.
%
%   ___ = CLIM(ax, ___ ) uses the axes specified by ax instead of the
%   current axes.
% 
%   CLIM sets or gets the CLim or CLimMode property of an axes.
%
%   See also COLORBAR, AXIS, YLIM, ZLIM, THETALIM, RLIM.

%   Copyright 2021 The MathWorks, Inc.


narginchk(0,2);
nargoutchk(0,2);

axHandle = [];
modeFlag = '';
climits = [];

if (nargin == 0)
    % Valid no arg syntaxes:
    %   CLIM
    %   vec = CLIM
    %   [min, max] = CLIM

    axHandle = gca;

elseif nargin == 1
    % Valid single arg syntaxes:
    %   CLIM('manual')
    %   CLIM('auto')
    %   M = CLIM('mode')
    %   CLIM(V)
    %   CLIM(AX)
    if (isempty(varargin{1}))
        error(message('MATLAB:caxis:InvalidVector'));
    elseif (matlab.graphics.internal.isCharOrString(varargin{1}))
        modeFlag = varargin{1};
        axHandle = gca;
    elseif (isscalar(varargin{1}) && isgraphics(varargin{1}))
        % Check for isscalar is necessary, else some limits values may be
        % interpreted as valid graphics handles, e.g [0 1]. 
        axHandle = varargin{1};
    else
        % otherwise, assume inputs are limits
        climits = varargin{1};
        axHandle = gca;
    end
elseif nargin == 2
    % Valid 2 arg syntaxes:
    %   CAXIS(AX, 'manual')
    %   CAXIS(AX, 'auto')
    %   M = CLIM(AX, 'mode')
    %   CAXIS(AX, V)
    if isgraphics(varargin{1})
        axHandle = varargin{1};
    else
        % All valid 2-arg syntaxes require initial axes handle.
        error(message('MATLAB:caxis:InvalidFirstArgument'))
    end

    if (isempty(varargin{2}))
        error(message('MATLAB:caxis:InvalidVector'))
    elseif(matlab.graphics.internal.isCharOrString(varargin{2}))
        modeFlag = varargin{2};
    else
        % otheriwse, assume inputs are limits
        climits = varargin{2};
    end
end

% Chart subclass support
if isa(axHandle,'matlab.graphics.chart.Chart')
    % If caxis is called with a chart handle, it will dispatched to the 
    % chart method directly; for this reason, varargin will never contain 
    % a chart handle. Here, we explictly provide the handle to invoke the
    % chart method in the gca case.
    try
        [varargout{1:nargout}] = clim(axHandle, varargin{:});
    catch me
        throw(me)
    end
    return
end

% Axes input must be scalar/non-empty & either cartesian, polar or 
% geographic axes.
if ~isscalar(axHandle) || ~( isgraphics(axHandle, 'axes') ...
        || isgraphics(axHandle, 'polaraxes') ...
        || isgraphics(axHandle, 'geoaxes'))
    error(message('MATLAB:caxis:NeedScalarAxesHandle'));
end

% Process cases where an output is populated.
if isempty(modeFlag) && isempty(climits)
    % Valid output-requested syntaxes:
    %   CAXIS     % ans populated in this case
    %   CAXIS(AX) % ans populated in this case
    %   vec = CAXIS
    %   vec = CAXIS(AX)
    %   [min, max] = CAXIS
    %   [min, max] = CAXIS(AX)
    currentCLims = get(axHandle,'CLim');
    if(nargout <= 1)
        varargout{1} = currentCLims;
    elseif(nargout == 2)
        varargout{1} = currentCLims(1);
        varargout{2} = currentCLims(2);
    end
    return
end

if ~isempty(modeFlag) 
    if strcmp(modeFlag, 'mode') 
        varargout{1} = get(axHandle, 'CLimMode');
    else
        nargoutchk(0,0);
        matlab.graphics.internal.markFigure(axHandle);
        try
            set(axHandle,'CLimMode',modeFlag);
        catch e
            throw(e);
        end
    end
elseif ~isempty(climits)
    if (numel(climits) == 2)
        set(axHandle, 'CLim',  climits);
    else
        error(message('MATLAB:caxis:InvalidNumberElements'))
    end
end
end
