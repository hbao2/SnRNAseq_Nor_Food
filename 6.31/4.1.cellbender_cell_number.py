import os
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import warnings
warnings.simplefilter("ignore", FutureWarning)
warnings.simplefilter("ignore", UserWarning)
warnings.simplefilter("ignore", RuntimeWarning)


# read cell bender output
file_dir = '/run/user/1000/gvfs/sftp:host=129.81.246.235/media/XLStorage/hbao2/SCRNA/Cellranger_output_20241020_6.31/clean_adata/'
adatas = [x for x in os.listdir(file_dir) if x.endswith('filtered.h5')]

def load_it(adata):
    samp = adata.split('_denoised')[0]
    adata = sc.read_10x_h5(file_dir + adata)
    adata.obs['Sample'] = samp
    adata.obs['dataset'] = 'TH_FCA'
    adata.obs.index = adata.obs.index + '_' + samp
    df = pd.read_csv('/home/denglab/Project/20241022_Host_Nor_Food/Gene_ID.txt', sep='\t', header=None, names=['FBgn', 'Symbol'])
    # Reorder df to match the order of adata.var_names
    df = df.set_index('FBgn').reindex(adata.var_names).reset_index() 
    # Set the 'FBgn' column as the index to match `adata.var_names
    adata.var['Symbol'] = df['Symbol'].values
    adata.var_names = adata.var['Symbol']    
    return adata

adatas = [load_it(ad) for ad in adatas]

# Extract cell counts and sample names
cell_counts = [adata.n_obs for adata in adatas]
sample_names = [adata.obs['Sample'][0] for adata in adatas]  # Assuming all rows in a sample have the same label

# Create a DataFrame for plotting
df = pd.DataFrame({'Sample': sample_names, 'Cell Count': cell_counts})

# Plot using seaborn
plt.figure(figsize=(8, 5))
sns.barplot(data=df, x='Sample', y='Cell Count', palette='viridis')

# Add labels and title
plt.title('Cell Numbers from Each Sample', fontsize=16)
plt.xlabel('Sample', fontsize=14)
plt.ylabel('Cell Count', fontsize=14)
plt.xticks(rotation=45)
plt.tight_layout()

# Show the plot
plt.show()
