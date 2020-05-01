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
