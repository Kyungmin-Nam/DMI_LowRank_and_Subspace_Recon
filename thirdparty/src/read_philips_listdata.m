function [data,info] = read_philips_listdata(name, Parameter2Read, scan_type)
% READ_PHILIPS_LIST_DATA     Reads a Philips .LIST/.DATA file pair
%
% [DATA,INFO] = READ_PHILIPS_LIST_DATA(NAME,PARAMETER2READ)
%
%   NAME is a string containing a file prefix or name of the .LIST
%   header file or .DATA binary file, e.g. RAW_001 or RAW_001.LIST or
%   RAW_001.DATA.
%
%   PARAMETER2READ is a structure containing the range of values used for
%   selectively extracting the data from the .DATA binary file. The 
%   PARAMETER2READ structure is filled using the information from the GI 
%   struture obtained by read_philips_list.m.
%
%   DATA is an N-dimensional array holding the raw binary data.
%
%   INFO is a structure containing details from the .LIST header file.
%
% Example:
%  [data,info] = readListData('raw_777.list');
%
%  After modification by Namgyun Lee
%  Parameter2Read.typ   = 'STD';
%  Parameter2Read.mix   = [];
%  Parameter2Read.dyn   = 0;
%  Parameter2Read.card  = 0;
%  Parameter2Read.echo  = 0;
%  Parameter2Read.loca  = (0:1).';
%  Parameter2Read.chan  = (0:1).';
%  Parameter2Read.extr1 = 0;
%  Parameter2Read.extr2 = 0;
%  Parameter2Read.aver  = 0;
%  [data,info] = read_philips_list_data('raw_777.list',Parameter2Read);
%
%  See also: LOADLABRAW, READ_PHILIPS_LIST
%
%  Dependencies: none
%

% Revision History
% * 2006.01.01    initial version - brianwelch
% * 2015.05.20    minor modification - Namgyun Lee
% * 2015.07.08    major modification and renamed as read_philips_listdata - Namgyun Lee
% modification to read spectro MRSI data by Namgyun Lee 2016-09-15

%------------------------------------------------------------------------
% Parse the option
%------------------------------------------------------------------------
if ~isempty(Parameter2Read.typ),   typ   = Parameter2Read.typ;   else typ   = 'STD'; end
if ~isempty(Parameter2Read.mix),   mix   = Parameter2Read.mix;   else mix   = 0;     end
if ~isempty(Parameter2Read.dyn),   dyn   = Parameter2Read.dyn;   else dyn   = 0;     end
if ~isempty(Parameter2Read.card),  card  = Parameter2Read.card;  else card  = 0;     end
if ~isempty(Parameter2Read.echo),  echo  = Parameter2Read.echo;  else echo  = 0;     end
if ~isempty(Parameter2Read.loca),  loca  = Parameter2Read.loca;  else loca  = 0;     end
if ~isempty(Parameter2Read.chan),  chan  = Parameter2Read.chan;  else chan  = 0;     end
if ~isempty(Parameter2Read.extr1), extr1 = Parameter2Read.extr1; else extr1 = 0;     end
if ~isempty(Parameter2Read.extr2), extr2 = Parameter2Read.extr2; else extr2 = 0;     end
if ~isempty(Parameter2Read.aver),  aver  = Parameter2Read.aver;  else aver  = 0;     end

toks = regexp(name, '^(.*?)(\.list|\.data)?$', 'tokens');
prefix = toks{1}{1};
list_name = sprintf('%s.list', prefix);
data_name = sprintf('%s.data', prefix);

%--------------------------------------------------------------------------
% Read a .LIST header file
%--------------------------------------------------------------------------
fid = fopen(list_name, 'r');
if fid ~= -1
    tstart = tic; fprintf('%s: Reading .LIST header file: %s ...', datetime, list_name);
    list_text = fread(fid, inf, 'uint8=>char').';
    fclose(fid);
    fprintf('done (%4.4f sec)\n', toc(tstart));
else
    error('cannot open %s for reading', list_name);
end

%--------------------------------------------------------------------------
% Set the INFO structure which contains details from the .LIST header file
%--------------------------------------------------------------------------
tstart = tic; fprintf('%s: Parsing the header file for %s data type ...', datetime, typ);
if strcmp(scan_type ,'spectro')
    attributes = {'mix' , 'dyn', 'card', 'echo', 'loca', 'chan', 'extr1', 'extr2', ...
                  'kx'  , 'ky' , 'kz'  , 'aver', 'sign', 'rf'  , 'grad' , 'enc'  , ...
                  'rtop', 'rr' , 'size', 'offset'};
else
    attributes = {'mix' ,  'dyn', 'card', 'echo', 'loca', 'chan', 'extr1', 'extr2', ...
                  'ky'  ,  'kz' , 'n.a.', 'aver', 'sign', 'rf'  , 'grad' , 'enc'  , ...
                  'rtop',  'rr' , 'size', 'offset'};
end
              
% Clean attributes names (will be used as variablenames and fieldnames)
for idx = 1:length(attributes)
    attributes{idx} = cleanFieldname(attributes{idx});
end

pattern = sprintf('(?<typ>%s+)', typ);
for idx = 1:length(attributes)
    range_flag = false;
    
    % Only for those attributes listed below a user will be able to select 
    % the subset of the total data
    switch attributes{idx}
        case 'mix'
            range_str = num2str(mix);   range_flag = true;
        case 'dyn'
            range_str = num2str(dyn);   range_flag = true;
        case 'card'
            range_str = num2str(card);  range_flag = true;
        case 'echo'
            range_str = num2str(echo);  range_flag = true;
        case 'loca'
            range_str = num2str(loca);  range_flag = true;
        case 'chan'
            range_str = num2str(chan);  range_flag = true;
        case 'extr1'
            range_str = num2str(extr1); range_flag = true;
        case 'extr2'
            range_str = num2str(extr2); range_flag = true;
        case 'aver'
            range_str = num2str(aver);  range_flag = true;
    end
    
    if range_flag
       range_str = num2str(range_str);
       range_str = range_str(~isspace(range_str)); % remove any blanks between characters

       pattern = sprintf('%s\\s+(?<%s>-?[%s]+)', pattern, attributes{idx}, range_str);
%        if strcmp(attributes{idx}, 'echo')
%           pattern = sprintf('%s\\s+(?<%s>-?%s)', pattern, attributes{idx}, range_str);
%        else
%           pattern = sprintf('%s\\s+(?<%s>-?[%s]+)', pattern, attributes{idx}, range_str);
%        end
    else
       pattern = sprintf('%s\\s+(?<%s>-?\\d+)', pattern, attributes{idx});
    end
end

info = regexp(list_text, pattern, 'names');
fprintf('done (%4.4f sec)\n', toc(tstart));

%--------------------------------------------------------------------------
% DATA is a multi-dimensional array organized as
% STD = Standard data vector (image data or spectroscopy data)
% REJ = Rejected standard data vector
%       (only for scans with arrhythmia rejection)
% NOI = Preparation phase data vector for noise determination
% NAV = Phase navigator data vector
% DNA = Dynamic phase navigator vector
% From 'correction .LIST/.DATA'
% PHX = Correction data vector for EPI/GraSE phase correction
% FRX = Correction data vector for frequency spectrum correction
% From 'measured .LIST/.DATA'
% PHC = Preparation phase data vector for calculation of phase correction
%       (only for EPI/GraSE scans), also known as 'echo-phase data'
% FRC = Preparation phase data vector for calculation of frequency response
%       correction
%--------------------------------------------------------------------------
% DATA is a multi-dimensional array organized as
if strcmp(scan_type, 'spectro')
    % kf is the readout (measurement) direction
    switch typ
        case 'STD'
            order = {'kx', 'ky', 'kz', 'kf', 'loca', 'dyn', 'card', 'echo', 'mix', 'aver', 'chan'};
        case 'REJ'
            order = {'kx', 'ky', 'kz', 'kf', 'loca', 'dyn', 'card', 'echo', 'mix', 'aver', 'chan'};
        case {'PHX','PHC'}
            order = {'kx', 'ky', 'kz', 'kf', 'loca', 'echo', 'mix', 'chan', 'sign', 'rf', 'grad'};
        case {'FRX','FRC'}
            order = {'kx', 'ky', 'kz', 'kf', 'loca', 'echo', 'mix', 'chan', 'sign'};
        case 'NOI'
            order = {'kf', 'loca', 'chan'};
        case 'NAV'
            order = {'kx', 'ky', 'kz', 'kf', 'loca', 'dyn', 'card', 'echo', 'mix', 'aver', 'chan'};
        case 'DNA'
            order = {'kf', 'dyn', 'chan'};
        otherwise
            order = {''};
    end
else
    % kx is readout (measurement) direction
    switch typ 
        case 'STD'
            order = {'kx', 'ky', 'kz', 'loca', 'dyn', 'card', 'echo', 'mix', 'aver', 'chan'};
        case 'REJ'
            order = {'kx', 'ky', 'kz', 'loca', 'dyn', 'card', 'echo', 'mix', 'aver', 'chan'};
        case {'PHX','PHC'}
            order = {'kx', 'ky', 'kz', 'loca', 'echo', 'mix', 'chan', 'sign', 'rf', 'grad'};
        case {'FRX','FRC'}
            order = {'kx', 'ky', 'kz', 'loca', 'echo', 'mix', 'chan', 'sign'};
        case 'NOI'
            order = {'kx', 'loca', 'chan'};
        case 'NAV'
            order = {'kx', 'ky', 'kz', 'loca', 'dyn', 'card', 'echo', 'mix', 'aver', 'chan'};
        case 'DNA'
            order = {'kx', 'dyn', 'chan'};
        otherwise
            order = {''};
    end
end

%--------------------------------------------------------------------------
% Determine the range of attributes of the selected type of complex data vector
%--------------------------------------------------------------------------
attributes_range.ky   = [];
attributes_range.kz   = [];
attributes_range.loca = [];
attributes_range.dyn  = [];
attributes_range.card = [];
attributes_range.echo = [];
attributes_range.mix  = [];
attributes_range.aver = [];
attributes_range.chan = [];
attributes_range.size = [];

% Determine the indices of the selected complex data type
idxTYPE = find(strcmp({info(:).typ},typ));

for idx = 1:length(order)
    if sum(strcmp(attributes,order{idx})) == 1
        list = str2num(char(unique({info(idxTYPE).(order{idx})})));
        attributes_range.(order{idx}) = [min(list) max(list)];
    end
end
attributes_range.size = str2num(char(unique({info(idxTYPE).size})));

%--------------------------------------------------------------------------
% Determine the size of a DATA matrix
%--------------------------------------------------------------------------
nr_attributes = length(order);
data_size = zeros(1,nr_attributes);

if strcmp(scan_type, 'spectro')
    redout_direction = 'kf';
else
    redout_direction = 'kx';
end

for idx = 1:nr_attributes
    if strcmp(order{idx}, redout_direction) == 1
        data_size(idx) = attributes_range.size/4/2;
    elseif strcmp(order{idx}, 'sign') == 1
        data_size(idx) = 2;
    else
        data_size(idx) = max(attributes_range.(order{idx})) - min(attributes_range.(order{idx})) + 1;
    end
end

% Preallocate data
data = zeros(data_size, 'single');

%--------------------------------------------------------------------------
% Open a .DATA binary file
%--------------------------------------------------------------------------
fid = fopen(data_name, 'r', 'ieee-le');
if fid == -1
    error('cannot open %s for reading', list_name);
end

hwait = waitbar(0,'=========================================================================================');
set(get(findobj(hwait, 'type', 'axes'), 'Title'), 'Interpreter', 'none');
set(get(findobj(hwait, 'type', 'axes'), 'Title'), 'String', sprintf('Reading %s raw data from %s ...', typ, data_name));
tstart = tic; fprintf('%s: Reading %s binary data: %s ...', datetime, typ, data_name);

N = length(idxTYPE);
for idx = 1:N
    if (fseek(fid, str2double(info(idxTYPE(idx)).offset), 'bof') == 0)
        % The complex data vectors are represented as binary data in little
        % endian single precision IEEE float format
        % (1 complex element = 2 floats = 8 bytes)
        tmp_data_1d = fread(fid, str2double(info(idxTYPE(idx)).size)/4, 'float32');
        tmp_data = tmp_data_1d(1:2:end) + 1j * tmp_data_1d(2:2:end);
        %tmp_data = tmp_data * str2double(info(idxTYPE(idx)).sign);

        tmpstr = '';
        for k = 1:length(order)
            switch order{k}
                case redout_direction % imaging: 'kx', spectro: 'kf'
                    tmpstr = sprintf('%s,1:%d', tmpstr, length(tmp_data));
                case 'sign'
                    if strcmp(info(idxTYPE(idx)).(order{k}),'-1')
                        tmpstr = sprintf('%s,1', tmpstr);
                    else
                        tmpstr = sprintf('%s,2', tmpstr);
                    end
                otherwise
                    index = str2double(info(idxTYPE(idx)).(order{k})) - ...
                            attributes_range.(order{k})(1) + 1;
                    tmpstr = sprintf('%s,%d', tmpstr, index);
            end
        end

        tmpstr(1) = []; % Delete initial comma
        eval(sprintf('data(%s) = tmp_data;', tmpstr));
    else
        error('Cannot FSEEK to offset=%d in data file %s', info(idxTYPE(k)).offset, data_name);
    end
    if mod(idx,100) == 99
       waitbar(idx/N, hwait);
    end
end

fclose(fid);
close(hwait);
fprintf('done (%4.4f sec)\n', toc(tstart));

end

function s = cleanFieldname(s)

illegal_chars = {...
    '+','-','*','.',...
    '^','\','/','.',...
    '=','~','<','>',...
    '&','|',':',';',...
    '(',')','{','}',...
    '[',']','{','}',...
    '''','%',' ','!', ...
    '@','#','$','`',...
    '?',',','"',...
    };

general_replacement_char = '_';
firstchar_replacement_char = 'x'; % cannot be an underscore

for k = 1:length(illegal_chars)
    s = strrep(s,illegal_chars{k},general_replacement_char);
end

% first character cannot be a number
firstchar_code = double(s(1));
if ( (firstchar_code >= double('0')) && (firstchar_code <= double('9')) )
    s(1) = firstchar_replacement_char;
end

% first character cannot be an underscore
if (s(1) == '_')
    s(1) = firstchar_replacement_char;
end

end