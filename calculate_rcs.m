function [] = calculate_rcs(directory, num_of_objects, frequency)



% %===============================================================================
% % Verify that no parallel pool is initialized by creating/deleting a dummy pool
% %===============================================================================
% if ~isempty(gcp('nocreate'))
%     poolobj = gcp('nocreate');
%     delete(poolobj);
% end
% 
% 
% %===============================================================================
% % Start parallel pool
% %===============================================================================
% % Get number of workers depending on job type. Start parallel pool via "local"
% % cluster object.
% %
% % obtain number of tasks from Slurm environment variables
% cluster = getenv('SLURM_CLUSTER_NAME');
% if strcmp(cluster, 'inter')
%     % interactive job
%     nw = str2num(getenv('SLURM_JOB_CPUS_PER_NODE'));
% elseif strcmp(cluster, 'cm2_tiny') || ...
%        strcmp(cluster, 'cm2') || ...
%        strcmp(cluster, 'serial') || ...
%        strcmp(cluster, 'mpp3')
%     % batch job
%     nw = str2num(getenv('SLURM_NTASKS_PER_NODE'));
% else
%     % default
%     nw = 28;
% end
% 
% % disallow Threading
% if maxNumCompThreads > 1
%     maxNumCompThreads(1);
%     warning('MultiThreading: number of threads has been set to 1!');
% end
% 
% % create a local cluster object
% pc = parcluster('local');
% % set number of workers
% pc.NumWorkers = nw;
% % set the JobStorageLocation to SCRATCH (default: HOME -> not recommended)
% pc.JobStorageLocation = strcat(getenv('SCRATCH'));
% % start the parallel pool
% poolobj = parpool(pc, nw);
% 
% 
% spmd
%     fprintf('Hello from MATLAB process PID=%d running on node %s!\n',...
%             feature('getpid'),...
%             getenv('HOSTNAME'));
% end


f = frequency;

save_path = "C:\Users\fatih\Desktop\new_objs\rcs\"

files = dir(directory);

stl_files = {};

for i = 1:length(files)
    % Check if the file is a .obj file
    [~, ~, ext] = fileparts(files(i).name);
    if strcmp(ext, '.stl')
        % Append the path of the .obj file to the list
        stl_files{end+1} = fullfile(files(i).folder, files(i).name);
    end
end


len = length(stl_files);
stl_files = string(stl_files);


files = stl_files(num_of_objects(1):len);
%files = stl_files(1:100);


az = -180:1:180;
for i = 1:numel(files)
    RCS_timetrial = zeros(361,181);
    filePath = files(i);
    p = platform;
    p.FileName = filePath;
    p.Units = "cm";
    
    [~,name,~] = fileparts(filePath);
    path = fullfile(save_path, sprintf('%s.mat',name));

    if isfile(path)
	   continue;
    else
        parfor c = 0:180
            el = c-90;
            X = sprintf('current elevation %d degrees',el);
            disp(filePath)
            disp(X)
            [values, ~, ~] = rcs(p,f,az,el,EnableGPU=0, Solver="PO", Polarization="HH");
            RCS_timetrial(:,c+1) = values;
        end
        R = RCS_timetrial;
        save(path,"RCS_timetrial","az")
    end

    %"/dss/dsshome1/03/ge26xih2/up_test/rcs_files"
    %fullsavepath = fullfile(save_path, name);
    %writematrix(a(i,:,:), fullsavepath)
end
end