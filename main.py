import os.path as op
import numpy as np      # array e operazioni numeriche veloci
import pandas as pd     # dataframe per enorme raccolta dati
import networkx as nx   # grafi efficienti e visualizzabili
import copy
import matplotlib.pyplot as plt # plotting

# output: dataframe con i dati dei pazienti
def get_PatientData(studyId):
    df = pd.DataFrame()
    path = f'data sets/{studyId}/data_clinical_patient.txt'
    if op.isfile(path):
        df = pd.read_csv(path, sep='\t', skiprows=4)
        df.drop_duplicates('PATIENT_ID', inplace=True)

    return df

# output: dataframe con i dati dei samples
def get_SampleData(studyId):
    df = pd.DataFrame()
    path = f'data sets/{studyId}/data_clinical_sample.txt'
    if op.isfile(path):
        df = pd.read_csv(path, sep='\t', skiprows=4)
        df.drop_duplicates(['SAMPLE_ID', 'PATIENT_ID'], inplace=True)

    return df

# output: dataframe con i dati delle mutazioni
def get_MutationData(studyId):
    df = pd.DataFrame()
    path = f'data sets/{studyId}/data_mutations.txt'
    if op.isfile(path):
        df = pd.read_csv(path, sep='\t')
        df.drop_duplicates(['Hugo_Symbol', 'Chromosome', 'Start_Position', 'End_Position', 'Tumor_Sample_Barcode'], inplace=True)

    return df

# output: dataframe con i dati dei trattamenti
def get_TreatmentData(studyId):
    df = pd.DataFrame()
    path = f'data sets{studyId}/data_timeline_treatment.txt'
    if op.isfile(path):
        df = pd.read_csv(path, sep='\t')

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
def build_PMGraph(full_data, patient_data, sample_data):
    graph = nx.from_pandas_edgelist(full_data, source='PATIENT_ID', target='MUTATION', create_using=nx.DiGraph())

    ps_data = pd.merge(patient_data, sample_data, left_on='PATIENT_ID', right_on='PATIENT_ID')
    ps_data.drop_duplicates('PATIENT_ID', inplace=True)
    attrs = ps_data.set_index('PATIENT_ID').to_dict('index')
    nx.set_node_attributes(graph, attrs)

    return graph

# output: lista dei nodi Di
def get_DiNodes(dip_graph):
    Di = [node for node, degree in dip_graph.in_degree() if degree == 0]

    return Di

# output: lista dei nodi P
def get_PNodes(pm_graph):
    pt = [(node, degree) for node, degree in pm_graph.out_degree() if degree != 0]
    P = [n for n, d in sorted(pt, key=lambda item: item[1], reverse=True)]

    return P

# output: lista dei nodi M
def get_MNodes(pm_graph):
    M = [node for node, degree in pm_graph.out_degree() if degree == 0]

    return M

# output: dataframe con ogni mutazione legata alla malattia in input con il numero di pazienti in cui compare
def getMutations_fromDisease(dip_graph, pm_graph, disease, mcount):
    col_name = 'Count'
    dcount = {}
    for p in dip_graph.neighbors(disease):
        for m in pm_graph.neighbors(p):
            if m in dcount:
                dcount[m] += 1
            else:
                dcount[m] = 1
    mutations_count = dict(sorted(dcount.items(), key=lambda x: x[1], reverse=True))
    mutations_df = pd.DataFrame(list(mutations_count.items()), columns=['Mutation', col_name])
    mutations_df['Frequency (%)'] = mutations_df.apply(lambda row: round((row[col_name] / mcount) * 100), axis=1)

    return mutations_df

# calculate similarity between two mutation sets
def cluster_similarity(mutations1, mutations2):
    common_mutations = mutations1 & mutations2
    all_mutations = mutations1 | mutations2
    s = len(common_mutations) / len(all_mutations)
    return s

def clustering(pm_graph, threshold=1):

    # clustering algorithm
    patients = get_PNodes(pm_graph)
    clusters = {}
    cc = 0
    for p in patients:
        mutations = set(m for m in pm_graph.neighbors(p))
        cluster_found = False

        for cl_number, cl_patients in clusters.items():
            cl_leader = cl_patients[0]
            leader_mutations = set(m for m in pm_graph.neighbors(cl_leader))
            similarity = cluster_similarity(mutations, leader_mutations)
            if similarity >= threshold:
                cl_patients.append(p)
                cluster_found = True
                break
        
        if not cluster_found:
            clusters[cc] = [p]
            cc += 1

    # sorting clusters and deleting those with just one patient
    clusters = dict(sorted(clusters.items(), key=lambda item: len(item[1]), reverse=True))
    final_clusters = {}
    cc = 0
    for n, patients in clusters.items():
        if len(patients) > 1:
            final_clusters[cc] = patients
            cc += 1
    
    return final_clusters

def main():
    return 0

if __name__ == "__main__":
    main()