import re
import subprocess

EXPANSIONS = {
    'acme1': {
        'machine_name': 'Acme1',
        'activation_path': '/usr/local/e3sm_unified/envs/load_latest_e3sm_unified_acme1.sh',
        'obs_path': '/p/user_pub/e3sm/e3sm_diags_data/obs_for_e3sm_diags/',
        'test_data_path': '/p/user_pub/e3sm/e3sm_diags_data/test_model_data_for_acme_diags/',
        'html_path': '/var/www/acme/acme-diags/<username>/',
        'web_address': 'https://acme-viewer.llnl.gov/<username>/',
    },
    'lcrc': {
        'machine_name': 'LCRC',
        'activation_path': '/lcrc/soft/climate/e3sm-unified/load_latest_e3sm_unified_chrysalis.sh',
        'obs_path': '/lcrc/soft/climate/e3sm_diags_data/obs_for_e3sm_diags/',
        'test_data_path': '/lcrc/soft/climate/e3sm_diags_data/test_model_data_for_acme_diags/',
        'html_path': '/lcrc/group/e3sm/public_html/diagnostic_output/<username>/',
        'web_address': 'https://web.lcrc.anl.gov/public/e3sm/diagnostic_output/<username>/',
    },
    'compy': {
        'machine_name': 'Compy',
        'activation_path': '/share/apps/E3SM/conda_envs/load_latest_e3sm_unified_compy.sh',
        'obs_path': '/compyfs/e3sm_diags_data/obs_for_e3sm_diags/',
        'test_data_path': '/compyfs/e3sm_diags_data/test_model_data_for_acme_diags/',
        'html_path': '/compyfs/www/<username>/',
        'web_address': 'https://compy-dtn.pnl.gov/<username>/',
    },
    'cori': {
        'machine_name': 'Cori',
        'activation_path': '/global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_cori-haswell.sh',
        'obs_path': '/global/cfs/cdirs/e3sm/e3sm_diags/obs_for_e3sm_diags/',
        'test_data_path': '/global/cfs/cdirs/e3sm/e3sm_diags/test_model_data_for_acme_diags/',
        'html_path': '/global/cfs/cdirs/e3sm/www/<username>/',
        'web_address': 'http://portal.nersc.gov/cfs/e3sm/<username>/'
    }
}


# Important: Do not have more than one `#expand` operation per line. The code will only expand the first one.
def generate_quick_guides():
    machine_names = EXPANSIONS.keys()
    git_top_level = subprocess.check_output(
        'git rev-parse --show-toplevel'.split()).strip().decode('utf-8')
    quick_guides_dir = '{g}/docs/source/quickguides/'.format(g=git_top_level)
    generic_quick_guide = '{d}quick-guide-generic.rst'.format(
        d=quick_guides_dir)
    specific_quick_guides = {}
    specific_quick_guide_files = {}
    for machine_name in machine_names:
        rst_file = '{d}quick-guide-{m}.rst'.format(
            d=quick_guides_dir, m=machine_name)
        specific_quick_guides[machine_name] = rst_file
        specific_quick_guide_files[machine_name] = open(rst_file, 'w')
    with open(generic_quick_guide, 'r') as file_read:
        # For line in file, if line matches #expand <expansion_name>#, then replace text with <expansion>.
        # Send output to the corresponding specific quick guide file, instead of overwriting the current file.
        for line in file_read:
            match_object = re.search('#expand ([^#]*)#', line)
            for machine_name in machine_names:
                multiprocessing = False
                if match_object is None:
                    new_line = line
                else:
                    expansion_indicator = match_object.group(0)
                    expansion_name = match_object.group(1)
                    if expansion_name == 'multiprocessing':
                        multiprocessing = True
                        new_line = ''
                    else:
                        expansion = EXPANSIONS[machine_name][expansion_name]
                        new_line = line.replace(expansion_indicator, expansion)
                if multiprocessing and machine_name != 'acme1':
                    multiprocessing_file = '{d}quick-guide-multiprocessing-{m}.rst'.format(
                        d=quick_guides_dir, m=machine_name)
                    with open(multiprocessing_file, 'r') as mpf:
                        for mpf_line in mpf:
                            specific_quick_guide_files[machine_name].write(
                                mpf_line)
                else:
                    specific_quick_guide_files[machine_name].write(new_line)
    for machine_name in machine_names:
        specific_quick_guide_files[machine_name].close()


if __name__ == '__main__':
    generate_quick_guides()
