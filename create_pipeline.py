import argparse
import configparser
import os


def write_batch(filename, script_path, path_to_write, script, script_args, container_path, *args):
    """
    write a single batch file from a provided python script
    :param script_path: path to python script
    :param path_to_write: path to write batch file
    :param filename: name of file to write to
    :param script: python script to be run in sbatch file
    :param script_args: any arguments the script requires as string with items separated by space
    :param container_path: path to container to use
    :param args: parameters for the cluster
    :return: sbatch file
    """
    params = {'script': script,
              'script_args': script_args,
              'script_path': script_path,
              'container_path': container_path,
              'n_tasks': args[0],
              'nodes': args[1],
              'job_name': args[2]}

    content = """#!/bin/bash\n
    #SBATCH --job-name={job_name}
    #SBATCH --nodes={nodes}
    #SBATCH --ntasks-per-node={n_tasks}
    #SBATCH --mem=4GB
    #SBATCH --time=00:10:00
    #SBATCH --output=logs/{job_name}-%j-stdout.log
    #SBATCH --error=logs/{job_name}-%j-stderr.log
    
    echo "Submitting SLURM job: {script} using {n_tasks} cores"
    mpirun singularity exec {container_path} python {script_path}{script} {script_args}
    """

    # if no script arguments provided, replace empty list with empty string
    if not params['script_args']:
        params['script_args'] = ''

    # insert arguments and remove whitespace
    content = content.format(**params).replace("    ", "")

    # if directory path doesn't exist create it
    if not os.path.exists(path_to_write):
        os.makedirs(path_to_write)
    write_location = path_to_write + '/' + filename
    sbatch_file = open(write_location, 'a')
    sbatch_file.write(content)
    sbatch_file.close()


def write_batch_all(scripts, script_path, script_args, container_path, *args):
    """
    write all the batch files
    :param script_path: python script path
    :param scripts: list of scripts
    :param script_args: list or arguments
    :param container_path: path to container to be used
    :param args: sbatch arguments
    :return: list of batch file names
    """
    script_list = []
    for task in args[0]:
        for i, script in enumerate(scripts):
            filename = str(script).replace('.py', '') + '_{}_{}_sbatch.sh'.format(str(task), i + 1)
            script_list.append('{}_cores'.format(str(task)) + '/' + filename)
            write_batch(filename, script_path, '{}_cores'.format(str(task)), script, script_args, container_path, task,
                        args[1], args[2])

    return script_list


def write_pipeline(script_list, filename):
    """
    write the final pipeline script submitted to slurm
    :param script_list: list of batch scripts (returned by write_batch_all() function)
    :param filename: name of file to be written
    :return:
    """
    pipeline_script = open(filename, 'a')
    pipeline_script.write('#!/bin/bash\n')
    pipeline_script.write('\n')
    pipeline_script.write('jid0=$(sbatch ' + script_list[0] + ')\n')
    pipeline_script.write('jid0=${jid0##* }\n')
    pipeline_script.write('$(echo "${jid0##* }" >> job_ids.txt)\n')
    for i, script in enumerate(script_list):
        if (i + 1) < len(script_list):
            pipeline_script.write('jid' + str(i + 1) + '=$(sbatch --dependency=afterok:$jid' + str(i) + ' ' +
                                  script_list[i + 1] + ')\n')
            pipeline_script.write('jid' + str(i + 1) + '=${jid' + str(i + 1) + '##* }\n')
            pipeline_script.write('$(echo "${jid' + str(i + 1) + '##* }" >> job_ids.txt)\n')

    pipeline_script.close()


def parse_configs(config_file, repetitions):
    """
    parse the config file for needed parameters
    :param config_file:
    :param repetitions: number of times to repeat single script if needed
    :return: tuple of parameters
    """
    config = configparser.ConfigParser()
    config.read(config_file)

    # parse sections of config
    scripts_information = config['script_information']
    slurm_information = config['slurm']

    # get script related information
    scripts = [eval(scripts_information['scripts'])] * repetitions
    script_path = eval(scripts_information['script_path'])
    script_args = eval(scripts_information['script_args'])
    container_path = eval(scripts_information['container_path'])

    # get slurm related information
    n_tasks = eval(slurm_information['n_tasks'])
    nodes = eval(slurm_information['nodes'])
    job_name = eval(slurm_information['job_name'])
    return scripts, script_path, script_args, container_path, n_tasks, nodes, job_name


def main(config_file, repetitions, pipeline_file):
    """
    run everything to create batch files and pipeline script
    :param config_file:
    :param repetitions:
    :param pipeline_file:
    :return:
    """
    py_scripts, py_script_path, py_script_args, sim_container_path, tasks_per_node, num_nodes, job_names \
        = parse_configs(config_file, repetitions)

    batch_list = write_batch_all(py_scripts, py_script_path, py_script_args, sim_container_path, tasks_per_node,
                                 num_nodes, job_names)
    write_pipeline(batch_list, pipeline_file)


if __name__ == '__main__':
    main('scripts_information.config', 1, 'pipeline_2.sh')
