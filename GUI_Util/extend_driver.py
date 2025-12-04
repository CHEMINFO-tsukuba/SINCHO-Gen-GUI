import sys, os, yaml
import subprocess
from glob import glob
import time
from IPython.core.debugger import Pdb

from ChemTSv2.generate_lead import Generate_Lead
from ChemTSv2.chemts_mothods import Methods, logs_dir
from AA_Score_Tool.MolEmbed import Embed_Mols
from AA_Score_Tool.ScoreCalc import AA_Score

import logging

if __name__ == '__main__':

    
    #### execution command checker ####
    args = sys.argv
    if len(args) != 2:
        print("Usage: python extend_driver.py ${yamlfile}")
        sys.exit(1)
    if not os.path.isfile(args[1]):
        print(f"{args[1]} does not exist!")
        sys.exit(1)
    #### execution command checker ####


    #### yaml loader ####
    config_file = str(sys.argv[1])
    with open(config_file, 'r')as f:
        config = yaml.safe_load(f)
    #####################

    
    #### working directory definition ####
    out_dir = config['OUTPUT']['directory']
    sampling_workflow_dir 	= os.path.join(out_dir, config['MD']['working_directory'])
    sincho_out_dir		= os.path.join(out_dir, config['SINCHO']['working_directory'])
    generate_working_directory 	= os.path.join(out_dir, config['ChemTS']['working_directory'])
    screening_dir 		= os.path.join(out_dir, config['AAScore']['working_directory'])
    logs_dir 			= os.path.join(out_dir, config['OUTPUT']['logs_dir'])

    os.makedirs(out_dir, exist_ok = True)
    os.makedirs(sampling_workflow_dir, exist_ok = True)
    os.makedirs(sincho_out_dir, exist_ok = True)
    os.makedirs(generate_working_directory, exist_ok = True)
    os.makedirs(screening_dir, exist_ok = True)
    os.makedirs(logs_dir, exist_ok = True)
    ######################################
    

    #### logger initialization ####
    logging.basicConfig(
        level 	= logging.INFO, 
        format 	= '%(asctime)s - %(levelname)s - %(message)s',
        stream	= sys.stdout
    )
    logger = logging.getLogger()
    file_handler = logging.FileHandler(os.path.join(logs_dir, 'main.log'))
    file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
    logger.addHandler(file_handler)
    
    logger.info(f"execution command 	: {args}")
    logger.info(f"output_dir  		: {out_dir}")
    logger.info(f"dir_tree 		: {sampling_workflow_dir}, {sincho_out_dir}, {generate_working_directory}, {screening_dir}, {logs_dir}")
    ###############################

    
    #### 01 Conformational Sampling using MD Simulation ####
    logger.info('01-Conformational Sampling Started')
    out_log_file = os.path.join(logs_dir, 'MD.log')
    with open(out_log_file, 'w') as stdout_f:
        subprocess.run(['python', '/MD-AUTOMATION/md_perform.py', config_file], stdout=stdout_f, stderr=stdout_f)
    ########################################################


    #### 02 Making Decision using SINCHO protocol ####
    logger.info('02-Making Decision Started')
    out_log_file = os.path.join(logs_dir, 'SINCHO.log')
    with open(out_log_file, 'w') as stdout_f:
        subprocess.run(['python', '/MD-AUTOMATION/p2c_sincho_parallel.py', config_file], stdout=stdout_f, stderr=stdout_f)
        subprocess.run(['python', '/SINCHO-H2L/yamlout.py', config_file], stdout=stdout_f, stderr=stdout_f)
    ##################################################
       
    trajectory_dirs = sorted(glob(os.path.join(sincho_out_dir, 'trajectory_*')))
 
    #### 03 Compound Generation using ChemTSv2 ####
    logger.info('ChemTS start.')
    generate_lead = Generate_Lead(trajectory_dirs, config)
    generate_lead.run()
    ###############################################
    

    ## AA Score
    ## AASCore only時の回避用
    if len(generate_lead.rank_output_dirs)==0:
        d_lists = []
        chemts_dir = os.path.join(config['OUTPUT']['directory'], config['ChemTS']['working_directory'])
        dirs = [ os.path.join(chemts_dir, f) for f in os.listdir(chemts_dir) if os.path.isdir(os.path.join(chemts_dir, f)) ]
        for traj_d in dirs:
            d_lists+=([ os.path.join(traj_d, d) for d in os.listdir(traj_d) if os.path.isdir(os.path.join(traj_d, d)) ])
        generate_lead.rank_output_dirs = d_lists
    print(generate_lead.rank_output_dirs)
    print(generate_lead.input_compound_files)


    logger.info('Mol Embed start.')
    embed_mols = Embed_Mols(generate_lead, config)
    embed_mols.run()

    logger.info('AA Score Calculation start.')
    aa_score = AA_Score(generate_lead, config)
    aa_score.run()

