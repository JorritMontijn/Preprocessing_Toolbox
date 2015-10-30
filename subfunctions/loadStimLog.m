function [sRec,intFlag] = loadStimLog(sRec)
	%UNTITLED2 Summary of this function goes here
	%   Detailed explanation goes here
	
	intFlag=0;
	strPath=['D:\Data\Raw\imaginglogs\' sRec.strSession filesep];
	try
		sRec.sStim = load([strPath sRec.sRawLib.strStimLog]);
		intFlag=1;
		sRec.sProcLog.boolStimLogFound = true;
	catch
		return;
	end
end

