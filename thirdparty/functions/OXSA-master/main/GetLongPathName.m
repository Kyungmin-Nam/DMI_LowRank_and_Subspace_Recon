function longFullFilename = GetLongPathName(varargin)
% Windows-only function that converts a supplied filename to the "long
% form" stored on disk, which includes normalising the case of the
% supplied filename.
%
% LongPathName = GetLongPathName(Name)
% INPUT:
%   Name: String or cell string, file or folder name with or without
%         relative or absolute path.
%         Unicode characters and UNC paths are supported.
%         Up to 1000 characters are allowed here, but some functions of the
%         operating system may support 260 characters only.
%
% OUTPUT:
%   FullName: String or cell string containing the normalised "long form"
%             file or folder name.
%
% EXAMPLES:
%   GetLongPathName({'C:\windows';'C:\'})
%
% $Id$

inFile = varargin{1};
%% Windows only

if ~ispc()
    warning([myFunc ':WindowsOnly'],'GetLongPathName is only effective on Microsoft Windows systems. Returning input unaltered.')
    if nargout > 0
        varargout{1} = varargin{1};
    end
    return
end


%% Load required libraries
myDir = fileparts(mfilename('fullpath'));

addpath(fullfile(myDir, 'getLongPathNameLibs'))

if libisloaded('kernel32'), unloadlibrary('kernel32'), end
%Requires C compiler:
% loadlibrary('kernel32', 'win_edited_ctr2.h','mfilename','loadkernel_ctr2')

%No compiler needed:
%#include \getLongPathNameLibs\kernel32_thunk_pcwin64.dll
%#include \getLongPathNameLibs\kernel32_thunk_pcwin64.exp
%#include \getLongPathNameLibs\kernel32_thunk_pcwin64.lib
%#include \getLongPathNameLibs\kernel32_thunk_pcwin64.obj
%#include \getLongPathNameLibs\win_edited_ctr2.h
loadlibrary('kernel32', @loadkernel_ctr2)

%Information on library functions
% libfunctions('kernel32')
% libfunctions('kernel32','-full')

%% Get full filename
%
unicode_input = typecast(unicode2native([inFile],'UTF-16'),'uint16');
unicode_input = [unicode_input(2:end) 0]; % remove BOM and NUL terminate
unicode_input2 = typecast(unicode2native(repmat(' ',1,1000),'UTF-16'),'uint16');
unicode_input2 = [unicode_input2(2:end) 0]; % remove BOM and NUL terminate


[a,b,c,d]=calllib('kernel32','GetFullPathNameW',...
    unicode_input,...
    1000,...
    unicode_input2,...
    []);


fullFilename = native2unicode(typecast([65279 c(1:a)],'uint8'),'UTF-16');

% N.B. This is still not guaranteed to have correct CASE...
% returnVal = lower(fullFilename)

%% GetLongPathName

unicode_input = typecast(unicode2native(fullFilename,'UTF-16'),'uint16');
unicode_input = [unicode_input(2:end) 0]; % remove BOM and NUL terminate
unicode_input2 = typecast(unicode2native(repmat(' ',1,1000),'UTF-16'),'uint16');
unicode_input2 = [unicode_input2(2:end) 0]; % remove BOM and NUL terminate
[a,b,c]=calllib('kernel32','GetLongPathNameW',...
    unicode_input,...
    unicode_input2,...
    1000);

if a ~= 0
    longFullFilename = native2unicode(typecast([65279 c(1:a)],'uint8'),'UTF-16');
else
    error('The specified directory does not exist.')
end

% N.B. This is still not guaranteed to have correct CASE...
longFullFilename = lower(longFullFilename);
