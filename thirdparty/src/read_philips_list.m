function [general_info, list_text] = read_philips_list(list_filename, scan_type)
% Written by Namgyun Lee
% Email: namgyunl@kbsi.re.kr, ggang56@gmail.com (preferred)
% Created: 07/08/2015, Last modified: 07/08/2015
% Updated by Namgyun Lee 05/04/2016 for MRSI
% Updated by Namgyun Lee 09/15/2016 for MRSI

%--------------------------------------------------------------------------
% Read a .LIST header file
%--------------------------------------------------------------------------
fid = fopen(list_filename,'r');
if fid ~= -1
    start_time = tic;
    fprintf('Reading .LIST header file: %s ...', list_filename);
    list_text = fread(fid,inf,'*char').';
    fclose(fid);
    elapsed_time = toc(start_time); 
    fprintf('done (%4.4f sec)\n', elapsed_time);
else
    error('cannot open %s for reading', list_filename);
end

%--------------------------------------------------------------------------
% Parse general information
%--------------------------------------------------------------------------
tokens_list = {
    % name of the field                  % number of values to be read
    %----------------------------------------------------------------------
    % number of ...
    %----------------------------------------------------------------------
    'number_of_mixes'                    , 1; ...
    %----------------------------------------------------------------------
    % number of ...
    %----------------------------------------------------------------------
    'number_of_encoding_dimensions'      , 1; ...
    'number_of_dynamic_scans'            , 1; ...
    'number_of_cardiac_phases'           , 1; ...
    'number_of_echoes'                   , 1; ...
    'number_of_locations'                , 1; ...
    'number_of_extra_attribute_1_values' , 1; ...
    'number_of_extra_attribute_2_values' , 1; ...
    'number_of_signal_averages'          , 1; ...
    %----------------------------------------------------------------------
    % number of ...
    %----------------------------------------------------------------------
    'number of coil channels'            , 1; ... % missing in MRS & MRSI
    %----------------------------------------------------------------------
    % k-space coordinate ranges
    %----------------------------------------------------------------------
    't_range'                            , 2; ... % MRS & MRSI
    'kx_range'                           , 2; ...
    'ky_range'                           , 2; ...
    'kz_range'                           , 2; ...
    %----------------------------------------------------------------------
    % k-space oversample factors
    %----------------------------------------------------------------------
    't_oversample_factor'                , 1; ... % MRS & MRSI
    'kx_oversample_factor'               , 1; ...
    'ky_oversample_factor'               , 1; ...
    'kz_oversample_factor'               , 1; ...
    %----------------------------------------------------------------------
    % number of ...
    %----------------------------------------------------------------------
    'number_of_rf_echoes'                , 1; ... % only for TSE, TFE, GraSE
    'number_of_gradient_echoes'          , 1; ... % only for EPI/GraSE
    %----------------------------------------------------------------------
    % reconstruction matrix
    %----------------------------------------------------------------------
    'F-resolution'                       , 1; ... % MRS & MRSI
    'X-resolution'                       , 1; ...
    'Y-resolution'                       , 1; ...
    'Z-resolution'                       , 1; ...
    %----------------------------------------------------------------------
    % SENSE factors (spatial dirs only!)
    %----------------------------------------------------------------------
    'X-direction SENSE factor'           , 1; ...
    'Y-direction SENSE factor'           , 1; ...
    'Z-direction SENSE factor'           , 1; ...
    %----------------------------------------------------------------------
    % imaging space coordinate ranges
    %----------------------------------------------------------------------
    'F_range'                            , 2; ... % MRS & MRSI
    'X_range'                            , 2; ...
    'Y_range'                            , 2; ...
    'Z_range'                            , 2; ...
    %----------------------------------------------------------------------
    % linear phase errors
    %----------------------------------------------------------------------
    '0th_order_phase_error_F'            , 1; ... % MRS & MRSI
    '1st_order_phase_error_F'            , 1; ... % MRS & MRSI
    '0th_order_phase_error_X'            , 1; ...
    '1st_order_phase_error_X'            , 1; ...
};

nr_filenames = length(tokens_list);

for idx = 1:nr_filenames
    nr_values = tokens_list{idx,2};
    if nr_values == 2
        pattern = sprintf('%s\\s+:\\s+([-0-9.]+\\s+[-0-9.]+)', tokens_list{idx});
    else
        pattern = sprintf('%s\\s+:\\s+([-0-9.]+)', tokens_list{idx});
    end
    tokens = regexp(list_text, pattern, 'tokens');

    nr_tokens = length(tokens);
    values = zeros(nr_tokens,nr_values);
    for idx2 = 1:nr_tokens
        values(idx2,:) = str2num(tokens{idx2}{:});
    end
    fieldname = cleanFieldname(tokens_list{idx});
    general_info.(fieldname) = values;
end

%--------------------------------------------------------------------------
% Parse the coil channel information (at the trailer of a .LIST file)
% NOTE: this is missing in the MRS & MRSI .LIST/.DATA file
%--------------------------------------------------------------------------
if ~strcmp(scan_type, 'spectro')
    fieldname = 'coil_channel_combination';
    pattern = sprintf('%s\\s+:\\s+([-0-9.]+)', fieldname);
    tokens = regexp(list_text, pattern, 'tokens');
    general_info.(fieldname) = tokens;

%     coil_channel_numbers = cell(general_info.number_of_locations,1);
%     for idx = 1:general_info.number_of_locations
%         coil_channel_numbers{idx} = strfind(general_info.coil_channel_combination{idx}{:},'1') - 1;
%     end
end

%--------------------------------------------------------------------------
% Remove empty fields in the general information
%--------------------------------------------------------------------------
field_names = fieldnames(general_info);
for idx = 1:nr_filenames
    if isempty(general_info.(field_names{idx}))
        general_info = rmfield(general_info, field_names{idx});
    end
end

end

function s = cleanFieldname(s)

illegal_chars = {'+',  '-', '*', '.', ...
                 '^',  '\', '/', '.', ...
                 '=',  '~', '<', '>', ...
                 '&',  '|', ':', ';', ...
                 '(',  ')', '{', '}', ...
                 '[',  ']', '{', '}', ...
                 '''', '%', ' ', '!', ...
                 '@',  '#', '$', '`', ...
                 '?',  ',', '"', ...
                 };

general_replacement_char = '_';
firstchar_replacement_char = 'x'; % cannot be an underscore

for idx = 1:length(illegal_chars)
    s = strrep(s,illegal_chars{idx},general_replacement_char);
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
