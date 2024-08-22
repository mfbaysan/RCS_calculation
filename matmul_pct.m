function [Cglob, comptime] = matmul_pct(size_A, size_B)
 
%===============================================================================
% MATLAB EXAMPLE: PARALLEL HELLO WORLD USING PCT TOOLBOX
%                 -> matrix-matrix multiplication C = A*B
%
% INPUT
%   size_A, size_B ... 2-element row vectors defining sizes of A and B
% OUTPUT
%   C ................ result
%   comptime ......... computation time (matrix product only)
%===============================================================================
 
%===============================================================================
% Check input
%===============================================================================
if nargin~=2
    error('Invalid number of input arguments!');
end
if size_A(2)~=size_B(1)
    error(sprintf('Dimension mismatch of A (%d columns) and B (%d rows)!',...
                  size_A(2), size_B(1)));
end
 
%===============================================================================
% Verify that no parallel pool is initialized by creating/deleting a dummy pool
%===============================================================================
if ~isempty(gcp('nocreate'))
    poolobj = gcp('nocreate');
    delete(poolobj);
end
 
%===============================================================================
% Start parallel pool
%===============================================================================
% Get number of workers depending on job type. Start parallel pool via "local"
% cluster object.
%
% obtain number of tasks from Slurm environment variables
cluster = getenv('SLURM_CLUSTER_NAME');
if strcmp(cluster, 'inter')
    % interactive job
    nw = str2num(getenv('SLURM_JOB_CPUS_PER_NODE'));
elseif strcmp(cluster, 'cm2_tiny') || ...
       strcmp(cluster, 'cm2') || ...
       strcmp(cluster, 'serial') || ...
       strcmp(cluster, 'mpp3')
    % batch job
    nw = str2num(getenv('SLURM_NTASKS_PER_NODE'));
else
    % default
    nw = 1;
end
 
% disallow Threading
if maxNumCompThreads > 1
    maxNumCompThreads(1);
    warning('MultiThreading: number of threads has been set to 1!');
end
 
% create a local cluster object
pc = parcluster('local');
% set number of workers
pc.NumWorkers = nw;
% set the JobStorageLocation to SCRATCH (default: HOME -> not recommended)
pc.JobStorageLocation = strcat(getenv('SCRATCH'));
% start the parallel pool
poolobj = parpool(pc, nw);
 
%===============================================================================
% Work
%===============================================================================
spmd
    fprintf('Hello from MATLAB process PID=%d running on node %s!\n',...
            feature('getpid'),...
            getenv('HOSTNAME'));
end
 
% generate well-defined matrices
NA = prod(size_A);
NB = prod(size_B);
A = reshape( linspace( 1,NA, NA), size_A );
B = reshape( linspace(NB, 1, NB), size_B );
 
% distribute data to workers and do parallel computation
spmd
    Aloc = codistributed(A, codistributor2dbc([nw 1]));
    Bloc = codistributed(B, codistributor2dbc([1 nw]));
    tic;
    Cloc = Aloc*Bloc;
    t = toc;
end
 
comptime = max(cell2mat(t(:)));
fprintf('parallel computation (MPI) of matrix-matrix product:\n');
fprintf('\tnumber of tasks (MATLAB workers) = %d\n', nw);
fprintf('\ttime = %.2f s\n', comptime);
 
Cglob = gather(Cloc);
 
%===============================================================================
% Close parallel pool
%===============================================================================
delete(poolobj);