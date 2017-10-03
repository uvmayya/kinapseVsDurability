function datacell_longpath = remove_shortpaths(datacell,minpathlength)
% return datacell with paths whose length is less than (or equal to) minpathlength removed 

holder = [];
for i = 1:length(datacell)
	if size(datacell{i},1)>minpathlength
		holder = [holder,i];
	end
end
datacell_longpath = datacell(holder);
