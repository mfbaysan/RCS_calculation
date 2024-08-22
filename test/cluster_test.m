myfunction = 'simple_func';
cluster_name = 'cm2_tiny';
partition_name = 'cm2_tiny';
walltime = '00:30:00';
tasks_per_node = 1;
num_worker = 16;
[job,ch] = job_config(myfunction, cluster_name, partition_name, walltime, tasks_per_node, num_worker);

% home dir : /dss/dsshome1/03/ge26xih2