obj = stlFileChecker("car4_scaled1.stl")
%scale object by at least 8 times
f = 94e9;
plat = platform(FileName="car4_scaled1.stl", Units="cm") % m, ft, in, cm, 
RCS_timetrial = zeros(361,181);
%object = readObj("sample3output.obj")
%tf = isWatertight(object)
%% 
figure
show(plat)
%obj = fegeometry("orion_nofbc.stl");
%figure
%show(p)
az = -180:1:180;

% for c = 0:180
%     el = c-90; % -90:1:90;
%     %el_c = el(c + 1); % Get the c'th element of el
%     X = sprintf('current elevation %d degrees',el);
%     disp(X)
%     [values, ~, ~] = rcs(plat,f,az,el, Solver="PO", Polarization="HH");
%     RCS_timetrial(:,c+1) = values;
% end
% 
% save("new_sample_rcs_plane.mat","RCS_timetrial","az");

% el = 0;
% tic
% rcs(plat,f,az,el, Solver="PO", Polarization="HH" ); %Solver="FMM", [values, azimuth, elevation] = 
% toc

directory = 'C:\Users\fatih\Desktop\new_objs\*\*';
interval_start=1;
interval_end=100;
frequency=94e9;

calculate_rcs(directory, [interval_start interval_end], frequency)
