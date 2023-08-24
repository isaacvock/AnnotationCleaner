'''
Automatically generate a mikado list
'''

import os

# Sample dictionary
scallop = {
    'use_scallop': True,
    'path': 'path/to/scallop/annotation',
    'label': 'sc',
    'strandedness': True,
    'score': 0,
    'is_reference': False,
    'exclude_redundant': False,
    'strip_cds': False
}


stringtie = {
    'use_stringtie': True,
    'path': 'path/to/stringtie/annotation',
    'label': 'st',
    'strandedness': True,
    'score': 0,
    'is_reference': False,
    'exclude_redundant': False,
    'strip_cds': False
}

provided_annotations = {
    'cl': {
        'path': 'path/to/class2/annotation',
        'strandedness': True,
        'score': 0,
        'is_reference': False,
        'exclude_redundant': False,
        'strip_cds': False
    },
    'sq':{
        'path': 'path/to/class2/annotation',
        'strandedness': True,
        'score': 0,
        'is_reference': False,
        'exclude_redundant': False,
        'strip_cds': False
    }
}


os.chdir("C:\\Users\\isaac\\Documents\\")

with open('output.txt', 'w') as f:

    if scallop["use_scallop"]:

        scallop.pop("use_scallop")

        row = [str(value) for value in scallop.values()]

        f.write('\t'.join(row) + '\n')
    
    if stringtie["use_stringtie"]:

        stringtie.pop("use_stringtie")

        row = [str(value) for value in stringtie.values()]

        f.write('\t'.join(row) + '\n')

    if ~bool(provided_annotations):
        
        for key, inner_dict in provided_annotations.items():

            row = [str(value) for value in inner_dict.values()]

            row.insert(1, key)

            f.write('\t'.join(row) + '\n')


    