import time # to get program's execution time
import pandas as pd
import numpy as np
from bravado.client import SwaggerClient # to access the API

# variable used in API functions for access the data
cb = SwaggerClient.from_url('https://www.cbioportal.org/api/v2/api-docs', config={"validate_requests":False,"validate_responses":False,"validate_swagger_spec": False})

# returns the list of all the cancer types (885)
def getAllCancerTypes():
    clist = []
    all_cancer_types = cb.Cancer_Types.getAllCancerTypesUsingGET().result()
    for c in all_cancer_types:
        clist.append(c.name)
    
    return clist

# given a study id and a disease, returns the list of all the samples for patients with that disease in the specified study
def getSamples(study_id, disease):
    samples = []
    all_samples = cb.Samples.getAllSamplesInStudyUsingGET(studyId = study_id).result()
    for s in all_samples:
        sample_id = s.sampleId
        data = cb.Clinical_Data.getAllClinicalDataOfSampleInStudyUsingGET(studyId = study_id, sampleId = sample_id, attributeId = 'CANCER_TYPE_DETAILED').result()
        if (data[0].value == disease):
            samples.append(sample_id)
    
    return samples

# returns the list of all the mutations for the specified gene
def getMutations(geneId, molecularProfile_id, sampleList_id, detail):
    mutations = cb.Mutations.getMutationsInMolecularProfileBySampleListIdUsingGET(
                    entrezGeneId = geneId,
                    molecularProfileId = molecularProfile_id,
                    sampleListId = sampleList_id,
                    projection = detail
                    ).result()
    
    return mutations

# build the first set of edges (diseases -> patients)
def build_dpGraph(associations, max_mc = None):
    #graph = {}
    df = pd.DataFrame(columns=['Paziente', 'Malattia'])
    for k, v in associations.items():
        disease = k
        study_id = v[0]
        samples = cb.Samples.getAllSamplesInStudyUsingGET(studyId = study_id).result()
        for s in samples:
            sample_id = s.sampleId
            data_d = cb.Clinical_Data.getAllClinicalDataOfSampleInStudyUsingGET(studyId = study_id, sampleId = sample_id, attributeId = 'CANCER_TYPE_DETAILED').result()
            d = data_d[0].value
            if (max_mc):
                data_mc = cb.Clinical_Data.getAllClinicalDataOfSampleInStudyUsingGET(studyId = study_id, sampleId = sample_id, attributeId = 'MUTATION_COUNT').result()
                if not data_mc:
                    continue
                mc = int(data_mc[0].value)
                if (mc > max_mc):
                    continue
            if (d == disease):
                p = data_d[0].patientId
                df.loc[len(df)] = [p, d]
                #if (d not in graph):
                    #graph[d] = [p]
                #elif (p not in graph[d]):
                    #graph[d].append(p)

    df.to_excel('dis_pat.xlsx', index=False)
    #return graph

# build the second set of edges (patients -> mutations)
def build_pmGraph(associations):
    #graph = {}
    dp_df = pd.read_excel('dis_pat.xlsx')
    genes = np.array([])
    with open("GreenGraph/genes.txt") as file:
        for line in file:
            x = line.strip()
            genes = np.append(genes, x)
    df = pd.DataFrame(columns=['Paziente', 'Malattia'])
    for g in genes:
        gene = cb.Genes.getGeneUsingGET(geneId = g).result()
        gene_id = gene.entrezGeneId
        for k, v in associations.items():
            molecularProfile_id = v[1]
            sampleList_id = v[2]
            mutations = getMutations(gene_id, molecularProfile_id, sampleList_id, "ID")
            for m in mutations:
                p = m.patientId
                if p in dp_df['Paziente'].values:
                    if p not in df['Paziente'].values:
                        df.loc[len(df), 'Paziente'] = p
                        df.loc[df['Paziente'] == p, 'Malattia'] = k
                    df.loc[df['Paziente'] == p, g] = 1
                #if p in dpGraph[k]:
                    #if (p not in graph):
                        #graph[p] = [g]
                    #elif (g not in graph[p]):
                        #graph[p].append(g)
    
    df.to_excel('cb_data.xlsx', index=False)
    #return graph

# get disgenet data for the analysed diseases
def getDisgenetData():
    data = []
    df = pd.read_excel('disgenet_data.xlsx')
    l = df.to_dict('records')
    for d in l:
        dt = dict(d)
        for k, v in d.items():
            if type(v) == float and np.isnan(v):
                del dt[k]
        data.append(dt)

    data1 = {}
    for d in data:
        for k, v in d.items():
            if k == 'Disease':
                disease = v
                data1[disease] = set()
            else:
                data1[disease].add(k)
    
    return data1

# get cbioportal data in dictionaries to work with it
def getCbData():
    data = []
    df = pd.read_excel('cb_data.xlsx')
    l = df.to_dict('records')
    for d in l:
        dt = dict(d)
        for k, v in d.items():
            if type(v) == float and np.isnan(v):
                del dt[k]
        data.append(dt)
    
    dpGraph = {}
    pmGraph = {}
    for d in data:
        for k, v in d.items():
            if (k == 'Paziente'):
                patient = v
                pmGraph[patient] = list()
            elif (k == 'Malattia'):
                disease = v
                if (disease not in dpGraph):
                    dpGraph[disease] = [patient]
                else:
                    dpGraph[disease].append(patient)
            else:
                pmGraph[patient].append(k)

    return (dpGraph, pmGraph)

def cb_toExcel(dpGraph, pmGraph):
    df = pd.DataFrame(columns=['Paziente', 'Malattia', 'Mutazioni'])
    for disease, patients in dpGraph.items():
        for patient in patients:
            genes = pmGraph.get(patient)
            if genes:
                df.loc[len(df)] = [patient, disease, genes]

    df.to_excel('cb_data (better view).xlsx', index=False)

def cb_common_toExcel(dpGraph, pmGraph):
    df = pd.DataFrame(columns=['Malattia', 'Mutazioni totali', 'Mutazioni in comune'])
    for disease in list(dpGraph):
        u = list(getMutationsfromDisease_union(dpGraph, pmGraph, disease))
        i = list(getMutationsfromDisease_intersect(dpGraph, pmGraph, disease))
        df.loc[len(df)] = [disease, u, i]

    df.to_excel('cb_data_common.xlsx', index=False)

def getMutationsfromDisease_union(dpGraph, pmGraph, disease):
    mutations = set()
    for p in dpGraph[disease]:
        if p in pmGraph:
            for m in pmGraph[p]:
                mutations.add(m)
    
    return mutations

def getMutationsfromDisease_intersect(dpGraph, pmGraph, disease):
    common_mutations = set()
    for p in dpGraph[disease]:
        if p in pmGraph:
            if not common_mutations:
                common_mutations = set(pmGraph[p])
            else:
                mutations = set(pmGraph[p])
                common_mutations = common_mutations.intersection(mutations)
    
    return common_mutations

def getMutationsfromDisease_percentage(dpGraph, pmGraph, disease):
    dp = {}
    for p in dpGraph[disease]:
        if p in pmGraph:
            for m in pmGraph[p]:
                if (m not in dp):
                    dp[m] = 1
                else:
                    dp[m] += 1

    n = len(dpGraph[disease])
    percentages = {}
    for m, v in dp.items():
        x = (v / n) * 100
        percentages[m] = round(x)
    sorted_d = sorted(percentages.items(), key=lambda x: x[1], reverse=True)
    
    return sorted_d

def getRelation(cb_data, ds_data):
    result = {}
    for k, v in cb_data.items():
        Md = ds_data[k]
        Mu = v[0]
        Mi = v[1]
        if (Md == Mi):
            result[k] = "Md = Mi : la conoscenza medica coincide perfettamente con l'evidenza medica"
        elif (Md.issubset(Mi)):
            result[k] = "Md < Mi : la conoscenza medica é incompleta ( ci sono dei geni mutati da ogni paziente con una determinata malattia che non sono noti essere legati a tale malattia"
        elif (Mi.issubset(Md)):
            result[k] = "Mi < Md : la conoscenza medica coincide parzialmente con l'evidenza medica ( ci sono dei geni legati a una determinata malattia che non sono mutati in tutti i pazienti con tale malattia)"
        else:
            result[k] = "Md e Mi sono scorrelati"
    return result

def allData_toExcel(dpGraph, pmGraph, ds_data):
    df2 = pd.DataFrame(columns=['Malattia', 'N.pazienti', 'Geni mutati da almeno un paziente', 'N. geni mutati da almeno un paziente', 
                                'Geni mutati da tutti i pazienti', 'Geni noti in Disgenet', 'N. geni noti in DisGeNET', 
                                'Geni mutati in comune tra i db', 'Confronto'])
    cb_data = {}
    for k in dpGraph:
        cb_data[k] = (getMutationsfromDisease_union(dpGraph, pmGraph, k), getMutationsfromDisease_intersect(dpGraph, pmGraph, k))
    relation = getRelation(cb_data, ds_data)
    for k, v in dpGraph.items():
        disease = k
        npatients = len(v)
        cb_mutations = cb_data[k][0]
        cb_nmutations = len(cb_data[k][0])
        common_mutations = cb_data[k][1]
        ds_mutations = ds_data[k]
        ds_nmutations = len(ds_data[k])
        common = cb_mutations.intersection(ds_mutations)
        result = relation[k]
        df2.loc[len(df2)] = [disease, npatients, cb_mutations, cb_nmutations, common_mutations, ds_mutations, ds_nmutations, common, result]
    
    df2.to_excel('all_data.xlsx', index=False)

def genesData_toExcel(dpGraph, pmGraph, ds_data, min_percentage = None):
    df = pd.DataFrame(columns=['Gene'])
    for disease in dpGraph:
        percentages = getMutationsfromDisease_percentage(dpGraph, pmGraph, disease)
        for m, v in percentages:
            if (min_percentage):
                if v < min_percentage:
                    continue
            if m in ds_data[disease]:
                flag = 'X'
            else:
                flag = 'V'
            if m not in df['Gene'].values:
                df.loc[len(df), 'Gene'] = m
            df.loc[df['Gene'] == m, (disease + ' (%)')] = (str(v) + ' - ' + flag)

    df.to_excel('genes_data.xlsx', index=False)

def main():
    start = time.time()
    associations = {'Uterine Carcinosarcoma/Uterine Malignant Mixed Mullerian Tumor' : ['ucs_tcga_pan_can_atlas_2018', 'ucs_tcga_pan_can_atlas_2018_mutations', 'ucs_tcga_pan_can_atlas_2018_all']}
                    #'Cervical Squamous Cell Carcinoma' : ['pancan_pcawg_2020', 'pancan_pcawg_2020_mutations', 'pancan_pcawg_2020_all'],
                    #'Pancreatic Adenocarcinoma' : ['paad_tcga_pan_can_atlas_2018', 'paad_tcga_pan_can_atlas_2018_mutations', 'paad_tcga_pan_can_atlas_2018_all'], 
                    #'Hepatocellular Carcinoma' : ['hcc_mskimpact_2018', 'hcc_mskimpact_2018_mutations', 'hcc_mskimpact_2018_all'], 
                    #'Gallbladder Cancer' : ['gbc_msk_2018', 'gbc_msk_2018_mutations', 'gbc_msk_2018_all']}
    
    #build_dpGraph(associations)
    #build_pmGraph(associations)
    ds_data = getDisgenetData()
    dpGraph, pmGraph = getCbData()
    #cb_toExcel(dpGraph, pmGraph)
    #cb_common_toExcel(dpGraph, pmGraph)
    #allData_toExcel(dpGraph, pmGraph, ds_data)
    genesData_toExcel(dpGraph, pmGraph, ds_data, 10)
    print("----------", round(time.time() - start, 2), "seconds ----------") # print the execution time

if __name__ == "__main__":
    main()