import os.path as op
import pandas as pd     # dataframe per raccolta e analisi dati
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
        df = pd.read_csv(path, sep='\t')
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

# output: dataframe dei dati sui farmaci
def get_DrugsData():
    df = pd.DataFrame()
    path = 'geni_farmaci.xls'
    if op.isfile(path):
        df = pd.read_excel(path)
        df.dropna(subset='Farmaci si/no', inplace=True)
        df.reset_index(drop=True, inplace=True)
    
    return df

# output: dataframe con i dati completi (paziente, malattia, mutazioni)
def get_FullData(sample_data, mutation_data):
    if (sample_data.empty or mutation_data.empty):
        full_data = pd.DataFrame()
    else:
        full_data = pd.merge(sample_data, mutation_data, left_on='SAMPLE_ID', right_on='Tumor_Sample_Barcode')
        full_data = full_data.astype(str)
        full_data['MUTATION'] = full_data[['Hugo_Symbol', 'Chromosome', 'Start_Position', 'End_Position']].agg('_'.join, axis='columns')
        #full_data['VAF'] = full_data['t_alt_count'] / (full_data['t_alt_count'] + full_data['t_ref_count']) * 100
    
    return full_data

# output: lista degli archi Di -> P (malattie -> pazienti)
def build_DiPGraph(full_data):
    graph = nx.from_pandas_edgelist(full_data, source='CANCER_TYPE_DETAILED', target='PATIENT_ID', create_using=nx.DiGraph())

    return graph

# output: lista degli archi P -> M (pazienti -> mutazioni)
def build_PMGraph(full_data, patient_data, sample_data):
    graph = nx.from_pandas_edgelist(full_data, source='PATIENT_ID', target='MUTATION', create_using=nx.DiGraph())

    ps_data = pd.merge(patient_data, sample_data, left_on='PATIENT_ID', right_on='PATIENT_ID')
    ps_data.drop_duplicates('PATIENT_ID', inplace=True)
    ps_data.astype(str)
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

#output: dataframe con ogni malattia e rispettivo numero di pazienti affetti
def show_diseases(dip_graph):
    diseases = get_DiNodes(dip_graph)
    df0 = {'DISEASE' : [], 'PATIENT COUNT' : []}
    for di in diseases:
        pcnt = len(set(p for p in dip_graph.neighbors(di)))
        df0['DISEASE'].append(di)
        df0['PATIENT COUNT'].append(pcnt)
    df = pd.DataFrame(df0)
    df.sort_values(by='PATIENT COUNT', ascending=False, inplace=True)
    df.reset_index(drop=True, inplace=True)

    return df

def view_disease(dip_graph, pm_graph, disease):
    df0 = {'PATIENT_ID' : [], 'MUTATION COUNT': [], 'SURVIVAL RATE': [], 'SURVIVAL STATUS': []}
    for p in dip_graph.neighbors(disease):
        df0['PATIENT_ID'].append(p)
        mut_cnt = len(set(m for m in pm_graph.neighbors(p)))
        df0['MUTATION COUNT'].append(mut_cnt)
        df0['SURVIVAL RATE'].append(pm_graph.nodes[p]['OS_MONTHS'])
        df0['SURVIVAL STATUS'].append(pm_graph.nodes[p]['OS_STATUS'])
    df = pd.DataFrame(df0)

    return df

# output: dataframe con ogni mutazione legata alla malattia in input con il numero di pazienti in cui compare
def getMutations_fromDisease(dip_graph, pm_graph, disease):
    pcnt = dip_graph.degree(disease)
    col_name = 'Count'
    dc = {}
    for p in dip_graph.neighbors(disease):
        for m in pm_graph.neighbors(p):
            if m in dc:
                dc[m] += 1
            else:
                dc[m] = 1
    dc_sorted = dict(sorted(dc.items(), key=lambda x: x[1], reverse=True))
    df = pd.DataFrame(list(dc_sorted.items()), columns=['Mutation', col_name])
    df['Frequency (%)'] = df.apply(lambda row: round((row[col_name] / pcnt) * 100, 1), axis=1)
    df['Gene'] = df.apply(lambda row: row['Mutation'].split('_')[0], axis=1)
    df = df[['Gene', 'Mutation', 'Count', 'Frequency (%)']]

    return df

def getGenes_fromDisease(dip_graph, pm_graph, disease):
    pcnt = dip_graph.degree(disease)
    col_name = 'Count'
    dc = {}
    for p in dip_graph.neighbors(disease):
        pgenes = set()
        for m in pm_graph.neighbors(p):
            gene = m.split('_')[0]
            pgenes.add(gene)
        for gene in pgenes:
            if gene in dc:
                dc[gene] += 1
            else:
                dc[gene] = 1
    df = pd.DataFrame(list(dc.items()), columns=['Gene', col_name])
    df['Frequency (%)'] = df.apply(lambda row: round((row[col_name] / pcnt) * 100, 1), axis=1)

    return df

def getGeneMutations_fromDisease(dip_graph, pm_graph, disease):
    pcnt = dip_graph.degree(disease)
    col_name = 'Count'
    dc = {}
    for p in dip_graph.neighbors(disease):
        for m in pm_graph.neighbors(p):
            gene = m.split('_')[0]
            if gene in dc:
                dc[gene] += 1
            else:
                dc[gene] = 1

# calcola la similaritá tra due insiemi di mutazioni
def cluster_similarity(pm_graph, patient1, patient2):
    p1_mutations = set(m for m in pm_graph.neighbors(patient1))
    p2_mutations = set(m for m in pm_graph.neighbors(patient2))
    common_mutations = p1_mutations & p2_mutations
    all_mutations = p1_mutations | p2_mutations
    s = len(common_mutations) / len(all_mutations)
    return s

# algoritmo per clusterizzare i pazienti in base alle mutazioni comuni
def clustering(pm_graph, threshold=1):

    # algoritmo di clustering
    patients = get_PNodes(pm_graph)
    clusters = {}
    cc = 0
    for p in patients:
        cluster_found = False
        for cl_number, cl_patients in clusters.items():
            cluster_found = True
            for clp in cl_patients:
                similarity = cluster_similarity(pm_graph, p, clp)
                if similarity < threshold:
                    cluster_found = False
                    break
            if cluster_found:
                cl_patients.append(p)
                break
    
        if not cluster_found:
            clusters[cc] = [p]
            cc += 1

    # ordina i clusters ed elimina quelli formati da un solo paziente
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