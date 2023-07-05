__author__ = "Johannes Köster"
__copyright__ = "Copyright 2016, Johannes Köster"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"


import tempfile
from pathlib import Path
from snakemake.shell import shell

samtools_opts = get_samtools_opts(snakemake)
extra = snakemake.params.get("extra", "")
shell(
    "scallop "
)