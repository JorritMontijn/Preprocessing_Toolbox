function structStim = getStimFramesPrePro(sRec)
	%UNTITLED2 Summary of this function goes here
	%   Detailed explanation goes here
	
	strField = fieldnames(sRec.sStim);
	if isfield(sRec.sStim.(strField{1}),'structStim')
		structStim = sRec.sStim.(strField{1}).structStim;
		if isfield(structStim,'ActOnPulses'), structStim.FrameOn = structStim.ActOnPulses;end
		if isfield(structStim,'ActOffPulses'), structStim.FrameOff = structStim.ActOffPulses;end
	elseif isfield(sRec.sStim.(strField{1}),'TrialNumber') && isfield(sRec.sStim.(strField{1}),'ActOnPulses') && isfield(sRec.sStim.(strField{1}),'ActOffPulses') %type 3
		structStim(:).Orientation = sRec.sStim.(strField{1}).Orientation;
		structStim(:).TrialNumber = sRec.sStim.(strField{1}).TrialNumber;
		structStim(:).FrameOn = sRec.sStim.(strField{1}).ActOnPulses;
		structStim(:).FrameOff = sRec.sStim.(strField{1}).ActOffPulses;
		if isfield(sRec.sStim.(strField{1}),'Contrast'), structStim(:).Contrast = sRec.sStim.(strField{1}).Contrast;end
		if isfield(sRec.sStim.(strField{1}),'SpatialFrequency'), structStim(:).SpatialFrequency = sRec.sStim.(strField{1}).SpatialFrequency;end
	elseif isfield(sRec.sStim.(strField{1}),'vecPresStimOri') && isfield(sRec.sStim.(strField{1}),'vecStimActOnFrames') %type 4
		%stim props
		structStim(:).Orientation = sRec.sStim.(strField{1}).vecPresStimOri;
		structStim(:).Contrast = sRec.sStim.(strField{1}).vecPresStimContrast;
		structStim(:).TrialNumber = 1:length(sRec.sStim.(strField{1}).vecTrialActStartFrames);
		
		%frame timing
		structStim(:).FrameTrialStart = sRec.sStim.(strField{1}).vecTrialActStartFrames;
		structStim(:).FrameOn = sRec.sStim.(strField{1}).vecStimActOnFrames;
		structStim(:).FrameOff = sRec.sStim.(strField{1}).vecStimActOffFrames;
		
		%secs timing
		structStim(:).SecsTrialStart = sRec.sStim.(strField{1}).vecTrialActStartSecs;
		structStim(:).SecsOn = sRec.sStim.(strField{1}).vecStimActOnSecs;
		structStim(:).SecsOff = sRec.sStim.(strField{1}).vecStimActOffSecs;
		
		%behavior
		structStim(:).cellResponse = sRec.sStim.(strField{1}).cellResponse;
		structStim(:).cellRespSecs = sRec.sStim.(strField{1}).cellRespSecs;
		structStim(:).cellRespFrames = sRec.sStim.(strField{1}).cellRespFrames;
		
		structStim(:).vecTrialResponse = sRec.sStim.(strField{1}).vecTrialResponse;
		structStim(:).vecTrialRespSecs = sRec.sStim.(strField{1}).vecTrialRespSecs;
		structStim(:).vecTrialRespFrames = sRec.sStim.(strField{1}).vecTrialRespFrames;
		
	elseif isfield(sRec.sStim.(strField{1}),'vecPresStimOri') %type 1
		structStim(:).Orientation = sRec.sStim.(strField{1}).vecPresStimOri;
		structStim(:).TrialNumber = 1:length(sRec.sStim.(strField{1}).vecStimActOnPulses);
		structStim(:).FrameOn = sRec.sStim.(strField{1}).vecStimActOnPulses;
		structStim(:).FrameOff = sRec.sStim.(strField{1}).vecStimActOffPulses;
	elseif isfield(sRec.sStim.(strField{1}),'matPresStimOri') %type 2
		structStim(:).Orientation = sRec.sStim.(strField{1}).matPresStimOri;
		structStim(:).TrialNumber = 1:length(sRec.sStim.(strField{1}).vecStimActOnPulses);
		structStim(:).FrameOn = sRec.sStim.(strField{1}).vecStimActOnPulses;
		structStim(:).FrameOff = sRec.sStim.(strField{1}).vecStimActOffPulses;
	else
		error([mfilename ':StimulationDataSchemeNotSupported'],'Stimulation file data was not recognized')
	end
end

