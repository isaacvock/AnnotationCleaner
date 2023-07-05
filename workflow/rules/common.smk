import glob

SAMP_NAMES = list(config['samples'].keys())

def get_input_bams(wildcards):
    return config["samples"][wildcards.sample]