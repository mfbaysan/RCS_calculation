
f = 95e9;

files = dir('D:\ShapeNet\shapenet-watertight\Bench_100\*\*');

stl_files = {};

for i = 1:length(files)
    % Check if the file is a .obj file
    [~, ~, ext] = fileparts(files(i).name);
    if strcmp(ext, '.stl')
        % Append the path of the .obj file to the list
        stl_files{end+1} = fullfile(files(i).folder, files(i).name);
    end
end

stl_files = string(stl_files);

%files = stl_files(num_of_objects[1]:num_of_objects[2]);
files = stl_files(14:15);


az = -180:1:180;
for i = 1:numel(files)
    RCS_timetrial = zeros(361,181);
    filePath = files(i);
    p = platform;
    p.FileName = filePath;
    p.Units = "m";
	parfor c = 0:180
		el = c-90;
		X = sprintf('current elevation %d degrees',el);
		disp(filePath)
		disp(X)
		[values, ~, ~] = rcs(p,f,az,el,EnableGPU=0);
		RCS_timetrial(:,c+1) = values;
	end
    R = RCS_timetrial;

    [~,name,~] = fileparts(filePath);
    

    %"/dss/dsshome1/03/ge26xih2/up_test/rcs_files"
    %fullsavepath = fullfile(save_path, name);
    %writematrix(a(i,:,:), fullsavepath)
end
