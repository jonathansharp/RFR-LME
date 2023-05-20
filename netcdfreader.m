function output = netcdfreader(path)

% extract variable names
info = ncinfo(path);
%extract data
for n = 1:size(info.Variables,2)
    % replace dashes with underscores
    temp_name = replace(info.Variables(n).Name,'-','_');
    output.(temp_name) = ncread(path,info.Variables(n).Name);
end

end