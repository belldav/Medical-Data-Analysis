import os.path as op
import numpy as np      # array e operazioni numeriche veloci
import pandas as pd     # dataframe per enorme raccolta dati
import networkx as nx   # grafi efficienti e visualizzabili
import matplotlib.pyplot as plt # plotting

# output: dataframe con i dati dei pazienti
def get_PatientData(studyId):
    df = pd.DataFrame()
    path = f'{studyId}/data_clinical_patient.txt'
    if op.isfile(path):
        df = pd.read_csv(path, sep='\t', skiprows=4)

    return df

# output: dataframe con i dati dei samples
def get_SampleData(studyId):
    df = pd.DataFrame()
    path = f'{studyId}/data_clinical_sample.txt'
    if op.isfile(path):
        df = pd.read_csv(path, sep='\t', skiprows=4)

    return df

# output: dataframe con i dati delle mutazioni
def get_MutationData(studyId):
    df = pd.DataFrame()
    path = f'{studyId}/data_mutations.txt'
    if op.isfile(path):
        df = pd.read_csv(f'{studyId}/data_mutations.txt', sep='\t', skiprows=2)
        df.drop_duplicates(['Hugo_Symbol', 'Chromosome', 'Start_Position', 'End_Position', 'Tumor_Sample_Barcode'], inplace=True)

    return df

# output: dataframe con i dati dei trattamenti
def get_TreatmentData(studyId):
    df = pd.DataFrame()
    path = f'{studyId}/data_timeline_treatment.txt'
    if op.isfile(path):
        df = pd.read_csv(f'{studyId}/data_timeline_treatment.txt', sep='\t')

    return df

# output: dataframe con i dati completi (paziente, malattia, mutazioni)
def get_FullData(sample_data, mutation_data):
    if (sample_data.empty or mutation_data.empty):
        full_data = pd.DataFrame()
    else:
        full_data = pd.merge(sample_data, mutation_data, left_on='SAMPLE_ID', right_on='Tumor_Sample_Barcode')
        full_data = full_data.astype(str)
        full_data['MUTATION'] = full_data[['Hugo_Symbol', 'Chromosome', 'Start_Position', 'End_Position']].agg('_'.join, axis='columns')
    
    return full_data

# output: lista degli archi Di -> P (malattie -> pazienti)
def build_DiPGraph(full_data):
    graph = nx.from_pandas_edgelist(full_data, source='CANCER_TYPE', target='PATIENT_ID', create_using=nx.DiGraph())

    return graph

# output: lista degli archi P -> M (pazienti -> mutazioni)
def build_PMGraph(full_data):
    graph = nx.from_pandas_edgelist(full_data, source='PATIENT_ID', target='MUTATION', create_using=nx.DiGraph())

    return graph

# output: lista dei nodi Di
def get_DiNodes(dip_graph):
    Di = [node for node, degree in dip_graph.in_degree() if degree == 0]

    return Di

# output: lista dei nodi P
def get_PNodes(pm_graph):
    P = [node for node, degree in pm_graph.in_degree() if degree == 0]

    return P

# output: lista dei nodi M
def get_MNodes(PM_graph):
    M = [node for node, degree in PM_graph.out_degree() if degree == 0]

    return M

# output: dataframe con ogni mutazione legata alla malattia in input con il numero di pazienti in cui compare
def getMutations_fromDisease(dip_graph, pm_graph, disease):
    npatients = dip_graph.degree(disease)
    col_name = f'Frequenza (su {npatients} pazienti)'
    dcount = {}
    for p in dip_graph.neighbors(disease):
        for m in pm_graph.neighbors(p):
            if m in dcount:
                dcount[m] += 1
            else:
                dcount[m] = 1
    mutations_count = dict(sorted(dcount.items(), key=lambda x: x[1], reverse=True))
    mutations_df = pd.DataFrame(list(mutations_count.items()), columns=['Mutazione', col_name])
    mutations_df['%'] = mutations_df.apply(lambda row: round((row[col_name] / npatients) * 100), axis=1)

    return mutations_df

def clustering(pm_graph, full_data):
    pnodes = get_PNodes(pm_graph)
    clusters = {}
    for p in pnodes:
        mutations = frozenset(n for n in pm_graph.neighbors(p))
        cluster_found = False

        for cl_mutations, cl_patients in clusters.items():
            if mutations == cl_mutations:
                cl_patients.add(p)
                cluster_found = True
                break
        
        if not cluster_found:
            clusters[mutations] = {p}

    clusters = dict(sorted(clusters.items(), key=lambda item: len(item[1]), reverse=True))
    cluster_dfs = {}
    cc = 0
    for k, v in clusters.items():
        if len(v) > 1:
            cl_data = {'Paziente': [], 'Malattia': []}
            mut_data = {'Mutazione' : []}
            
            for p in v:
                disease = full_data.loc[full_data['PATIENT_ID'] == p, 'CANCER_TYPE_DETAILED'].values[0]
                cl_data['Paziente'].append(p)
                cl_data['Malattia'].append(disease)
            cl_df = pd.DataFrame(cl_data)

            for m in k:
                mut_data['Mutazione'].append(m)
            mut_df = pd.DataFrame(mut_data)
            cluster_dfs[cc] = [cl_df, mut_df]
            cc += 1

    return cluster_dfs

def main():
    return 0

if __name__ == "__main__":
    main()