import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def plot_mutation_frequency(site_list, gene_size, density=True, gene_name='', mutation_type=''):
    """To plot the mutation frequency of a gene

    Args:
        site_list (list of int): list containing sites of mutation
        gene_size (int): the size of the gene of interest
        density (bool, optional): Plot as density for cumulative frequency graph. Defaults to True.
        gene_name (string, optional): Name of the gene. Defaults to ''.
        mutation_type (string, optional): Type of mutation. Defaults to ''.
    """
    values, base = np.histogram(site_list, bins=gene_size, density=density)
    cumulative = np.cumsum(values)
    plt.plot(base[:-1], cumulative, c='blue')
    x_ticks = np.arange(0, gene_size, 100)
    plt.xticks(x_ticks)
    plt.xlim([0,gene_size])
    plt.title(f'{gene_name} {mutation_type} Mutation Frequency')
    plt.xlabel('Amino Acid Position')
    plt.ylabel('Cumulative frequency')

def plot_cancer_type_freq(cancer, gene_name='', mutation_type=''):
    """Plot the mutation frequency by cancer type

    Args:
        cancer (pandas series): perform a value_counts on the pandas df on cancer type
        gene_name (string, optional): Name of the gene. Defaults to ''.
        mutation_type (string, optional): Type of mutation. Defaults to ''.
    """
    plt.figure(figsize=(15,10)) 
    cancer.plot.bar()
    plt.title(f'{gene_name} {mutation_type} Mutation Frequency')
    plt.xlabel('Types of Cancer')
    plt.ylabel('Frequency of mutation')
    plt.tight_layout()

class Plot_site:
    def __init__(self):
        self.site_dict = {'Brain': [1150, 300],
                          'Breast': [3800, 2300],
                          'Liver': [1000, 2700],
                          'Cervix': [3400, 4250],
                          'Lung': [800, 2000],
                          'Uterus': [3400, 4000]}

    def get_site_dict(self):
        return self.site_dict

    def update_site_dict(self,site_dict):
        self.site_dict = site_dict

    def map_site(self, df):
        if 'site' not in df.columns:
            print("Error! There is no column 'site'. Please ensure the \
                  dataframe contains 'site' and 'frequency'")
        df['coordinate'] = df['site'].map(self.site_dict)
        return df

    def plot_site(self, df, img_path, size=2000, cmap=plt.cm.PiYG, figsize=(15,10)):
        if 'coordinate' not in df.columns:
            df = self.map_site(df)

        c_scale = 1/df.frequency.max()
        img = plt.imread(img_path)

        plt.figure(figsize=figsize)
        plt.imshow(img)
        for index, row in df.iterrows():
            plt.scatter(row['coordinate'][0], row['coordinate'][1], c=cmap(np.array(row['frequency']*c_scale).reshape(1,)), s=row['frequency']*size)
            plt.text(row['coordinate'][0]+200, row['coordinate'][1]+50, row['site'], fontsize=15, backgroundcolor='black', color='white')
            plt.text(row['coordinate'][0]-100, row['coordinate'][1]+50,f'{round(row.frequency*100)}%', fontweight=800)
        plt.tick_params(
            bottom=False,              
            left=False,
            labelleft=False,
            labelbottom=False) 
