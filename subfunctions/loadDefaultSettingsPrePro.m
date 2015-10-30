function sPS = loadDefaultSettingsPrePro()
	sLoad = load('PreProDefaultValues.mat');
	sPS = sLoad.sPS;
	clear sLoad;
end