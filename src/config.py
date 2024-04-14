import os


#------------------------------------------------------------------------------#
#                                NOTE TO USERS                                 #
#------------------------------------------------------------------------------#
# The three fields defined below are user- or system-specific. An email        #
# address is required to submit jobs to the NEOS server, which is how the      #
# nonlinear problems in `scripts/pre-06-model-options.ipynb` and               #
# `scripts/ijoc-01-lpac-geometries.ipynb` are solved. The other two fields     #
# should, by default, point to the `data` and `results` directories in this    #
# repository                                                                   #
#                                                                              #
EMAIL = 'username@provider.ext'
DATA_DIR = '/path/to/2023.0125/data'
RESULTS_DIR = '/path/to/2023.0125/results'
#                                                                              #
# The other definitions in this file need not be modified.                     #
#------------------------------------------------------------------------------#


# original data
ORIG_DATA_DIR = os.path.join(DATA_DIR, 'original')
USA_STATES_ZIP = os.path.join(ORIG_DATA_DIR, 'geoBoundaries-USA-ADM1-all.zip')
MEX_STATES_ZIP = os.path.join(ORIG_DATA_DIR, 'geoBoundaries-MEX-ADM1-all.zip')

# preprocessed data
PREP_DATA_DIR = os.path.join(DATA_DIR, 'preprocessed')
HEAD_TO_TAIL_JSON = os.path.join(PREP_DATA_DIR, 'head_to_tail.json')
SPECS_ACTIVS_2000_YAML = os.path.join(PREP_DATA_DIR, 'specs-ACTIVSg2000.yaml')
SPECS_REDUCED_ACTIVS_2000_YAML = os.path.join(PREP_DATA_DIR, 'specs-reduced-ACTIVSg2000.yaml')
SPECS_FLOOD_YAMLS = {
    ('imelda',): os.path.join(PREP_DATA_DIR, 'specs-imelda.yaml'),
    ('imelda', 'ev'): os.path.join(PREP_DATA_DIR, 'specs-imelda-ev.yaml'),
    ('imelda', 'mv'): os.path.join(PREP_DATA_DIR, 'specs-imelda-mv.yaml'),
    ('harvey',): os.path.join(PREP_DATA_DIR, 'specs-harvey.yaml'),
    ('harvey', 'ev'): os.path.join(PREP_DATA_DIR, 'specs-harvey-ev.yaml'),
    ('harvey', 'mv'): os.path.join(PREP_DATA_DIR, 'specs-harvey-mv.yaml')
}
SPECS_PF_OPTION_YAMLS = {
    'dc': os.path.join(PREP_DATA_DIR, 'specs-dc-options.yaml'),
    'lpacc': os.path.join(PREP_DATA_DIR, 'specs-lpacc-options.yaml'),
    'lpacf': os.path.join(PREP_DATA_DIR, 'specs-lpacf-options.yaml'),
    'qpac': os.path.join(PREP_DATA_DIR, 'specs-qpac-options.yaml')
}
SPECS_BIGMS_YAMLS = {
    'dc': os.path.join(PREP_DATA_DIR, 'specs-dc-bigM.yaml'),
    'lpacc': os.path.join(PREP_DATA_DIR, 'specs-lpacc-bigM.yaml'),
    'lpacf': os.path.join(PREP_DATA_DIR, 'specs-lpacf-bigM.yaml'),
    'qpac': os.path.join(PREP_DATA_DIR, 'specs-qpac-bigM.yaml')
}
SPECS_COMPLETE_YAMLS = {
    ('imelda', 'dc'): os.path.join(PREP_DATA_DIR, 'complete-imelda-dc.yaml'),
    ('imelda', 'dc', 'ev'): os.path.join(PREP_DATA_DIR, 'complete-imelda-ev-dc.yaml'),
    ('imelda', 'dc', 'mv'): os.path.join(PREP_DATA_DIR, 'complete-imelda-mv-dc.yaml'),
    ('imelda', 'lpacc'): os.path.join(PREP_DATA_DIR, 'complete-imelda-lpacc.yaml'),
    ('imelda', 'lpacc', 'ev'): os.path.join(PREP_DATA_DIR, 'complete-imelda-ev-lpacc.yaml'),
    ('imelda', 'lpacc', 'mv'): os.path.join(PREP_DATA_DIR, 'complete-imelda-mv-lpacc.yaml'),
    ('imelda', 'lpacf'): os.path.join(PREP_DATA_DIR, 'complete-imelda-lpacf.yaml'),
    ('imelda', 'lpacf', 'ev'): os.path.join(PREP_DATA_DIR, 'complete-imelda-ev-lpacf.yaml'),
    ('imelda', 'lpacf', 'mv'): os.path.join(PREP_DATA_DIR, 'complete-imelda-mv-lpacf.yaml'),
    ('imelda', 'qpac'): os.path.join(PREP_DATA_DIR, 'complete-imelda-qpac.yaml'),
    ('imelda', 'qpac', 'ev'): os.path.join(PREP_DATA_DIR, 'complete-imelda-ev-qpac.yaml'),
    ('imelda', 'qpac', 'mv'): os.path.join(PREP_DATA_DIR, 'complete-imelda-mv-qpac.yaml'),
    ('harvey', 'dc'): os.path.join(PREP_DATA_DIR, 'complete-harvey-dc.yaml'),
    ('harvey', 'dc', 'ev'): os.path.join(PREP_DATA_DIR, 'complete-harvey-ev-dc.yaml'),
    ('harvey', 'dc', 'mv'): os.path.join(PREP_DATA_DIR, 'complete-harvey-mv-dc.yaml'),
    ('harvey', 'lpacc'): os.path.join(PREP_DATA_DIR, 'complete-harvey-lpacc.yaml'),
    ('harvey', 'lpacc', 'ev'): os.path.join(PREP_DATA_DIR, 'complete-harvey-ev-lpacc.yaml'),
    ('harvey', 'lpacc', 'mv'): os.path.join(PREP_DATA_DIR, 'complete-harvey-mv-lpacc.yaml'),
    ('harvey', 'lpacf'): os.path.join(PREP_DATA_DIR, 'complete-harvey-lpacf.yaml'),
    ('harvey', 'lpacf', 'ev'): os.path.join(PREP_DATA_DIR, 'complete-harvey-ev-lpacf.yaml'),
    ('harvey', 'lpacf', 'mv'): os.path.join(PREP_DATA_DIR, 'complete-harvey-mv-lpacf.yaml'),
    ('harvey', 'qpac'): os.path.join(PREP_DATA_DIR, 'complete-harvey-qpac.yaml'),
    ('harvey', 'qpac', 'ev'): os.path.join(PREP_DATA_DIR, 'complete-harvey-ev-qpac.yaml'),
    ('harvey', 'qpac', 'mv'): os.path.join(PREP_DATA_DIR, 'complete-harvey-mv-qpac.yaml'),
}


# results
HEURISTIC_RESULTS_DIR = os.path.join(RESULTS_DIR, 'heuristic')
WS_RESULTS_DIR = os.path.join(RESULTS_DIR, 'ws')
EV_RESULTS_DIR = os.path.join(RESULTS_DIR, 'ev')
EEV_RESULTS_DIR = os.path.join(RESULTS_DIR, 'eev')
SP_RESULTS_DIR = os.path.join(RESULTS_DIR, 'sp')
SP_NG_RESULTS_DIR = os.path.join(RESULTS_DIR, 'sp-ng')  # ng = no good
MV_RESULTS_DIR = os.path.join(RESULTS_DIR, 'mv')
MMV_RESULTS_DIR = os.path.join(RESULTS_DIR, 'mmv')
RO_RESULTS_DIR = os.path.join(RESULTS_DIR, 'ro')
RO_NG_RESULTS_DIR = os.path.join(RESULTS_DIR, 'ro-ng')  # ng = no good
SPSOL_ROMOD_RESULTS_DIR = os.path.join(RESULTS_DIR, 'spsol-romod')
ROSOL_SPMOD_RESULTS_DIR = os.path.join(RESULTS_DIR, 'rosol-spmod')
ACPF_RESULTS_DIR = os.path.join(RESULTS_DIR, 'acpf')
CONSOLIDATED_RESULTS_DIR = os.path.join(RESULTS_DIR, 'consolidated')
FIG_RESULTS_DIR = os.path.join(RESULTS_DIR, 'figures')
