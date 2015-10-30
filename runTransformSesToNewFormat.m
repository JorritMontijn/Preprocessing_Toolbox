%function runTransformSesToNewFormat%( input_args )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% set session to convert
strDir = 'D:\Data\Processed\imagingdata\20130313\xyt02\';
strSes = '20130313xyt02_ses.mat';

%% perform conversion
%load data
strNewSes = [strSes(1:(end-4)) '2' strSes((end-3):end)];
sLoad=load([strDir strSes]);
ses = sLoad.ses;

%assign new additional field to all neurons
intNeurons = numel(ses.neuron);
for intNeuron=1:intNeurons
	ses.neuron(intNeuron).type = 'neuron';
end

%convert interneuron structures to neurons
cellInterneurons{1} = 'PV';
cellInterneurons{2} = 'SOM';
cellInterneurons{3} = 'VIP';

for intInterneuron=1:length(cellInterneurons)
	strType = cellInterneurons{intInterneuron};
	
	if isfield(ses,strType)
		intCells = numel(ses.(strType));
		for intCell=1:intCells
			intNeuron = intNeuron + 1;
			structCell = ses.(strType)(intCell);
			
			cellFields = fieldnames(structCell);
			for intField=1:length(cellFields)
				ses.neuron(intNeuron).(cellFields{intField}) = structCell.(cellFields{intField});
			end
			ses.neuron(intNeuron).type = strType;
		end
	end
end

%save new session file
save([strDir strNewSes],'ses')

