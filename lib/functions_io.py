### Imports ###
import pandas as pd

### Functions ###
def load_file(
    filename: str,
) -> pd.DataFrame:
    
    extension = filename.split('.')[-1]

    if extension == 'tsv':
        # Expression files
        df = pd.read_csv(filename, sep='\t', index_col=0, header=None)
    elif extension == 'feather':
        raw_df = pd.read_feather(filename)
        if 'gene1' in raw_df.columns:
            # Correlation networks
            df = raw_df.set_index(['gene1', 'gene2'])
        elif 'tf' in raw_df.columns and 'gene' in raw_df.columns:
            # PANDA networks
            df = raw_df.set_index(['tf', 'gene'])
        elif 'gene' in raw_df.columns:
            # Indegree files
            df = raw_df.set_index('gene')
        elif 'tf' in raw_df.columns:
            # Outdegree files
            df = raw_df.set_index('tf')
        else:
            # Assumed to be coexpression matrices
            df = raw_df.set_index('0')
    else:
        raise ValueError(f'Invalid extension for: {filename}, '
            'expected tsv or feather')

    return df