function ses = doTransformSpikes(input)
	%UNTITLED Summary of this function goes here
	%   Detailed explanation goes here
	
	%check input type
	if ischar(input)
		boolSave = true;
		sLoad = load(input);
		ses = sLoad.ses;
		clear sLoad;
	elseif isstruct(input)
		boolSave = false;
		ses = input;
		clear input;
	end
	
	%transform spikes
	intNeurons=numel(ses.neuron);
	for intNeuron=1:intNeurons
		vecSpikes = zeros(size(ses.neuron(intNeuron).dFoF));
		vecSpikes(ses.neuron(intNeuron).apFrames) = ses.neuron(intNeuron).apSpikes;
		ses.neuron(intNeuron).vecSpikes = vecSpikes;
	end
	
	%save file if needed
	if boolSave
		save(input,'ses');
	end
end

