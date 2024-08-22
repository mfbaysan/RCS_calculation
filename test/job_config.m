function [job,ch] = job_config(example, ...
                               cluster_name, ...
                               partition_name, ...
                               walltime, ...
                               tasks_per_node, ...
                               num_worker)
 
%===============================================================================
% MATLAB MPS EXAMPLE
% -> configuration script to initialze MPS and submit job
%===============================================================================
% input:
% example .................... name of user function without extension ".m"
% cluster_name ............... name of cluster, e. g.: cm2 (refers to the HPC
%                              machine, e.g. CoolMUC-2)
% partition_name ............. name of queue/partition, e. g.: cm2_std
%                              Details on cluster/partition names:
%                              https://doku.lrz.de/x/AgaVAg
% walltime, tasks_per_node ... equivalent to Slurm parameters
% num_worker ................. number of MPS workers
% return values:
% job ........................ job handle
% ch ......................... cluster object handle
%===============================================================================
 
%===============================================================================
% Step 1: cluster configuration
%===============================================================================
configCluster(cluster_name, partition_name);
 
%===============================================================================
% Step 2: job configuration
%===============================================================================
ch = parcluster;
ch.AdditionalProperties.WallTime = walltime;
ch.AdditionalProperties.ProcsPerNode = tasks_per_node;
ch.AdditionalProperties.AdditionalSubmitArgs = '--cpus-per-task=1 --mem=55G';


jobdir = fullfile(getenv('SCRATCH'), 'MdcsDataLocation/coolmuc/', version('-release'));
if ~exist(jobdir)
    mkdir(jobdir);
end
ch.JobStorageLocation = jobdir;


ch.saveProfile;
 
%===============================================================================
% Step 3: job submission to Slurm
%===============================================================================
% Command:
%   job = ch.batch(@myfunction, n_arg_out, {arg_in_1, ..., arg_in_n}, 'Pool', np)
% Input:
%   @myfunc ...... reference to user-defined function/script myfunc.m
%   n_arg_out .... number of expected output arguments
%   {arg_in_#} ... list of input arguments of myfunction
%   'Pool', np ... key-value-pair: define size of pool (number of workers)
%
% Help via Matlab commands:
%   help batch
%   doc batch
 
fhandle = eval(sprintf('@%s', example));
job = ch.batch(fhandle, 1, {}, 'Pool', num_worker);