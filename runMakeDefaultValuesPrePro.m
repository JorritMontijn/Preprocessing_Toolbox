
%make default values
sPS = struct;
for intRec=1:99
	sPS.strRecording{intRec} = sprintf('xyt%02d',intRec);
end
sPS.boolDoXMLRead = true; %true
sPS.boolDoStimLogRead = true; %true
sPS.boolDoTestLoad =  true; %true
sPS.boolDoSmooth = false; %false
sPS.strSmoothMethod = 'gaussian'; %'gaussian'
sPS.intSmoothKernelSize = 10; %10
sPS.boolDoRegistration = true; %true
sPS.intSubPixelDepth = 100; %100
sPS.boolDoAverageSaved = true; %true
sPS.boolSupervised = false; %false
sPS.boolDoRemSaturated = true;
sPS.boolDoPhaseCorrect = true;

save('PreProDefaultValues.mat','sPS')