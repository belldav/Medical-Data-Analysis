import time # to get program's execution time
import pandas as pd
import numpy as np
from bravado.client import SwaggerClient # to access the API

# variables used in API functions for access the data
cb = SwaggerClient.from_url('https://www.cbioportal.org/api/v2/api-docs',
                            config={"validate_requests":False,"validate_responses":False,"validate_swagger_spec": False})
studyId  = 'pancan_pcawg_2020'
moleculareProfileId = 'pancan_pcawg_2020_mutations'
sampleListId = 'pancan_pcawg_2020_all'

patients = set() # to store the analysed patients (106)
cancer_types = set() # to store the analysed cancer types (8)

# returns the list of all the cancer types (885)
def getAllCancerTypes():
    clist = []
    all_cancer_types = cb.Cancer_Types.getAllCancerTypesUsingGET().result()
    for c in all_cancer_types:
        clist.append(c)
    
    return clist

# returns the list of all the samples for patients with tumor stage >= 4 (112)
def getSamples(id):
    samples = []
    stages = ['4', 'IV', 'IVA', 'pT4aN0M0', 'pT4aN1M0', 'pT4aN2bM0', 'T4aN2', 'T4N0M0', 'T4N1', 'T4N1bM1', 'T4N1M0', 'T4N1M1', 'T4N2M0', 'T4N2M1', 'T4N3M0', 'T4NXMX']
    all_samples = cb.Samples.getAllSamplesInStudyUsingGET(studyId = id).result()
    for s in all_samples:
        sid = s.sampleId
        data = cb.Clinical_Data.getAllClinicalDataOfSampleInStudyUsingGET(studyId = id, sampleId = sid).result()
        for d in data:
            if (d.clinicalAttributeId == 'STAGE' and d.value in stages):
                samples.append(sid)
                break
    
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

# build the first set of edges (diseases -> patients) of the green graph
def build_dpGraph(id):
    graph = {}
    with open("GreenGraph/samples.txt") as file:
        samples = [line.strip() for line in file]
    for s in samples:
        data = cb.Clinical_Data.getAllClinicalDataOfSampleInStudyUsingGET(studyId = id, sampleId = s).result()
        for d in data:
            if (d.clinicalAttributeId == 'CANCER_TYPE'):
                disease = d.value
                patient = d.patientId
                break
        if (disease not in graph):
            graph[disease] = [patient]
        elif (patient not in graph[disease]):
            graph[disease].append(patient)
        patients.add(patient)
        cancer_types.add(disease)
    
    return graph

# build the second set of edges (patients -> mutations) of the green graph
def build_pmGraph():
    graph = {}
    genes = np.array([])
    with open("genes.txt") as file:
        for line in file:
            x = line.strip()
            genes = np.append(genes, x)
    for g in genes:
        gene = cb.Genes.getGeneUsingGET(geneId = g).result()
        geneId = gene.entrezGeneId
        mutations = getMutations(geneId, moleculareProfileId, sampleListId, "ID")
        for m in mutations:
            pid = m.patientId
            if pid in patients:
                if (pid not in graph):
                    graph[pid] = [g]
                elif (g not in graph[pid]):
                    graph[pid].append(g)
    
    return graph

def getMutationsfromDisease_union(dpGraph, pmGraph, disease):
    mutations = set()
    for p in dpGraph[disease]:
        for m in pmGraph[p]:
            mutations.add(m)
    
    return mutations

def getMutationsfromDisease_intersect(dpGraph, pmGraph, disease):
    common_mutations = set()
    for p in dpGraph[disease]:
        if not common_mutations:
            common_mutations = set(pmGraph[p])
        else:
            mutations = set(pmGraph[p])
            common_mutations = common_mutations.intersection(mutations)
    
    return common_mutations


def main():
    start = time.time()
    #disease = 'Ovarian Cancer'
    #dpGraph = build_dpGraph(studyId)
    #pmGraph = build_pmGraph()
    #m_union = getMutationsfromDisease_union(dpGraph, pmGraph, disease)
    #m_intersect = getMutationsfromDisease_intersect(dpGraph, pmGraph, disease)
    #print(m_union)
    #print(m_intersect)
    dpGraph = {'Ovarian Cancer': ['DO46336', 'DO46360', 'DO46382', 'DO46394', 'DO46493', 'DO46560'], 
               'Renal Cell Carcinoma': ['DO47243'], 
               'Pancreatic Cancer': ['DO48539', 'DO49183', 'DO51180', 'DO51187', 'DO221544', 'DO52123', 'DO52129', 'DO52132', 'DO52136', 'DO52137', 'DO52151', 'DO52159', 'DO52161', 'DO52164', 'DO32976', 'DO35132'], 
               'Hepatobiliary Cancer': ['DO48717', 'DO48738', 'DO50778', 'DO50780', 'DO50785', 'DO50793', 'DO50806', 'DO50842', 'DO50844', 'DO217864', 'DO218491', 'DO218535', 'DO218695', 'DO23520', 'DO23525', 'DO23529', 'DO23530', 'DO45092', 'DO45131', 'DO45193', 'DO45209', 'DO45221', 'DO45223', 'DO45281', 'DO45297', 'DO45305', 'DO45307'], 
               'Esophagogastric Cancer': ['DO50316', 'DO50321', 'DO50440', 'DO50447', 'DO50409', 'DO50411', 'DO50330', 'DO50337', 'DO50342', 'DO50354', 'DO50383', 'DO50336', 'DO50436', 'DO224498', 'DO217814', 'DO218442', 'DO217817', 'DO218487', 'DO222288', 'DO218443', 'DO217822', 'DO218223'], 
               'Prostate Cancer': ['DO52516'], 
               'Head and Neck Cancer': ['DO217836', 'DO217950', 'DO218012', 'DO218075', 'DO218180', 'DO218211', 'DO218280', 'DO218554', 'DO218697', 'DO218773'], 
               'Melanoma': ['DO220864', 'DO220883', 'DO220878', 'DO220879', 'DO220874', 'DO220875', 'DO220894', 'DO220877', 'DO220898', 'DO220899', 'DO220901', 'DO220902', 'DO220903', 'DO220904', 'DO220905', 'DO220906', 'DO220907', 'DO220908', 'DO220909', 'DO220910', 'DO220911', 'DO220912', 'DO220913']}
    df = pd.DataFrame(columns = ['Patient', 'Disease'])
    df.loc[len(df)] = ['Davide', 'Cancro']
    print(df)
    print("----------", round(time.time() - start, 2), "seconds ----------") # print the execution time

if __name__ == "__main__":
    main()
