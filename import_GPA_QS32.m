function [mass_index] = import_GPA_QS32(path)
%% Import_GPA for Quadstar32 v2 04/02/2021
%% by  Vincent Thor�ton
% This function imports any set of GPA data recorded using Quadstar32
% The function is called with a line such as: import_GPA(raw_data_path)
%It returns an array. The 1st row contains the recorded masses, the 2nd
%contains the corresponding column number.
global experiment_name experiment_date experiment_time possible_mass_num possible_mass_str raw_time raw_time_min MS_pressure raw_data_headers raw_data col reactor_pressure;

% [raw_data_filename,raw_data_root_folder] = uigetfile('*.*');
% path = strcat(raw_data_root_folder,raw_data_filename);
%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following gauges and masses may be imported:
possible_QMS_gauges={'PKR360','PKR251','PKR261','PKR','QMS pressure','DEFAULT'};
possible_reaction_gauges={'CMR','CMR261','CMR361','reactor pressure'};
possible_mass_str={'14','16','17','18','19','20','21','22','23','24','28','30','32','34','36','40','44','46','48'};%must be synchronised with the following:
possible_mass_num=str2double(possible_mass_str);%must be synchronised with the previous.
raw_data=reconstruct_data_file(sprintf('%s',path));
%juice the datafile:
experiment_name=raw_data{1,4};
experiment_date=raw_data{2,2};%experiment_date=raw_data(find(strcmp(raw_data,'DATE :')),2);%
experiment_time=raw_data{2,4};
%blocks_number=str2double(raw_data{7,4});%for future dev
%read what masses and pressure were recorded:
for ii=0:9
eval(sprintf('row=find(strcmp(raw_data,''Datablock %g''));',ii))
if isempty(row)==0
eval(sprintf('datablock_positions(%g)=row;',ii+1))
end
end
cycles_position=find(strcmp(raw_data,'Cycle'));

if size(datablock_positions,2)==1
    datablock_length(1)=cycles_position-2-datablock_positions(1);
elseif size(datablock_positions,2)==2
        datablock_length(1)=datablock_positions(2)-1-datablock_positions(1);
    datablock_length(2)=cycles_position-2-datablock_positions(2);
else
for ii=1:size(datablock_positions)-1
    datablock_length(ii)=datablock_positions(ii+1)-datablock_positions(ii)-1;
end
datablock_length(ii+1)=cycles_position-2-datablock_positions(ii+1);
end
shift=0;
for ii=1:size(datablock_positions,2)
raw_data(cycles_position,5+shift:shift+4+datablock_length(ii))=raw_data(datablock_positions(ii)+1:datablock_positions(ii)+datablock_length(ii),2:2)';
shift=shift+datablock_length(ii);
end

%delete the header:
raw_data(1:cycles_position-1,:) = [];
% delete the 3 1st columns:
raw_data(:,1:3) = [];
raw_data_headers=raw_data(1,:);%extract the content of the columns
raw_data(1,:) = [];%remove it
raw_data(end,:) = [];%remove the last row by precaution (could be NaN)
raw_data=str2double(raw_data);%transform strings to numbers

col=find(strcmp(raw_data_headers(1,:),'RelTime[s]'));%column 1 should be time; just checking
    raw_time=raw_data(:,col);
    raw_time_min=raw_time/60;

for ii=1:size(possible_QMS_gauges,2)
c=find(strcmp(raw_data_headers(1,:),possible_QMS_gauges(ii)));
if isempty(c)==0
col=c;
MS_pressure=raw_data(:,col);
if strcmp(possible_QMS_gauges{ii},'DEFAULT')
  MS_pressure=exp(MS_pressure/1000*1.667-11.33 );
end
end
end

for ii=1:size(possible_reaction_gauges,2)
c=find(strcmp(raw_data_headers(1,:),possible_reaction_gauges(ii)));
if isempty(c)==0
col=c;
reactor_pressure=raw_data(:,col);
end
end
raw_data_headers=string(round(str2double(raw_data_headers)));%round the values of the masses (with mass calibration the masses are not equal to unit
% locate the masses of interest and create a specific array for each of them
mass_index=[];%list of the recorded mass and corresponding column
for i = 1:size(possible_mass_str,2)
col=find(strcmp(raw_data_headers(1,:),possible_mass_str{1,i}));
if strcmp(raw_data_headers(1,col),possible_mass_str{1,i}) == 1%check i
mass_index=[mass_index [possible_mass_num(1,i);col]];
end
end
end