from collections import defaultdict
import numpy as np
from collections.abc import Mapping
import csv
import normalize_nodes
import pandas as pd

#We only need this as a class so that I can dynamically load in the npz for the MONDOs
#we're currently querying.
class ARAX_NPZ(Mapping):
    def __init__(self):
        self.active_npz_idx = None
        self.active_npz = 1
    def __getitem__(self,npz_idx):
        if(self.active_npz_idx!=npz_idx):
            del self.active_npz
            self.active_npz_idx = npz_idx
            self.active_npz = np.load(f"/home/ubuntu/ARAX/MONDO_NPZ/mondos_to_normalized_drugs_{npz_idx}.npz")["arr_0"]
        return self.active_npz
    def __iter__(self):
        return iter(range(0,22))
    def __len__(self):
        return 22
    

def buildDrugList(with_unnormed=False):
    arax_drug_curie_to_idx = {}
    arax_drug_curie_list = []
    if(with_unnormed):
        with open("/home/ubuntu/ARAX/KGML-xDTD/normalized_drug_idxs.csv") as f:
            reader = csv.reader(f)
            for i,row in enumerate(reader):
                arax_drug_curie_list.append(row[0:3])
                arax_drug_curie_to_idx[row[1]] = i
    else:
        with open("/home/ubuntu/ARAX/KGML-xDTD/normalized_drug_idxs.csv") as f:
            reader = csv.reader(f)
            for i,row in enumerate(reader):
                arax_drug_curie_list.append(row[1:3])
                arax_drug_curie_to_idx[row[1]] = i

    return (arax_drug_curie_to_idx, arax_drug_curie_list)

def buildMondoList():
    arax_mondo_to_idx = {}
    arax_mondo_list = []
    with open("/home/ubuntu/ARAX/KGML-xDTD/mondos_from_arax.txt") as f:
        for i,line in enumerate(f):
            arax_mondo_list.append(line.strip())
            arax_mondo_to_idx[line.strip()] = i
    return (arax_mondo_to_idx,arax_mondo_list)

class ValidArax():
    def __init__(self):
        (arax_mondo_to_idx,arax_mondo_list) = buildMondoList()
        self.arax_mondo_set = set(arax_mondo_list)
        self.arax_mondo_to_idx = arax_mondo_to_idx
        (arax_drug_curie_to_idx, arax_drug_curie_list) = buildDrugList()
        arax_drug_curie_set = set([x[0] for x in arax_drug_curie_list])
        self.arax_drug_curie_set = arax_drug_curie_set
    def getDrugsInARAX(self,drug_curie_iter):
        hits = []
        for drug_curie in drug_curie_iter:
            if(drug_curie in self.arax_drug_curie_set): hits.append(drug_curie)
        return hits
    def getMONDOsInARAX(self,mondo_iter):
        hits = []
        for mondo in mondo_iter:
            if(mondo in self.arax_mondo_set): hits.append(mondo)
        return hits




def buildTPData():
    (arax_mondo_to_idx,arax_mondo_list) = buildMondoList()
    (arax_drug_curie_to_idx, arax_drug_curie_list) = buildDrugList()
    tp_pairs = set()
    import csv
    tp_pairs = set()
    with open("/home/ubuntu/ARAX/arax_training_data.csv") as f:
        dreader = csv.DictReader(f)
        for row in dreader:
            disease_curie = row['normalized_disease_id']
            drug_curie = row['normalized_drug_id']
            if(disease_curie in arax_mondo_to_idx and drug_curie in arax_drug_curie_to_idx):
                disease_idx = arax_mondo_to_idx[disease_curie]
                drug_idx = arax_drug_curie_to_idx[drug_curie]
                tp_pairs.add((disease_idx,drug_idx))
    return tp_pairs

tp_pairs = buildTPData()
#This class builds the first layer of the dictonary. So ARAX_DICT["MONDO:1234567"] builds this.
class ARAX_MONDO_Map(Mapping):
    def __init__(self):
        self.arax_npz = ARAX_NPZ()
        (arax_mondo_to_idx,arax_mondo_list) = buildMondoList()
        self.arax_mondo_list = arax_mondo_list
        self.arax_mondo_to_idx = arax_mondo_to_idx
        (arax_drug_curie_to_idx, arax_drug_curie_list) = buildDrugList()
        self.arax_drug_curie_list = arax_drug_curie_list
        self.arax_drug_curie_to_idx = arax_drug_curie_to_idx
#        self.arax_mondo_list = None
#        self.arax_mondo_to_idx = None
        #self.arax_mondo_to_idx,self.arax_mondo_list = buildMondoList
        
    def __getitem__(self, key):
        return ARAX_Drug_Map(key,self)
        #return "abc"
    def __iter__(self):
        return iter(self.arax_mondo_list)
    def __len__(self):
        return len(self.arax_mondo_to_idx)
    
    #This function exists to simply enable us to quickly query over the top n 
    # top hits for a given MONDO. It basically only exists so that we don't
    # need to reimplement the code to utilize the dynamic loader for the npz
    # files. Since we'll almost certainly be hitting the MONDOs in a distinct pattern.
    def getTopNHitsForMONDOFromArax(self,mondo,n=500,check_for_tps=True):
        drug_map = ARAX_Drug_Map(mondo,self)
        np_list = drug_map.arax_npz[drug_map.mondo_bin][drug_map.mondo_sub_idx]
        #https://stackoverflow.com/questions/6910641/how-do-i-get-indices-of-n-maximum-values-in-a-numpy-array
        top_hit_indicies = np_list.argsort()[-n*2:][::-1]
        #ind = np.argpartition(np_list, -n)[-n:]
        results = []
        disease_idx = self.arax_mondo_to_idx[mondo]
        for drug_idx in top_hit_indicies:
            #This is 
            if(check_for_tps):
                if((disease_idx,drug_idx) in tp_pairs):continue
            #Tuples with (drug_curie, drug_name, score)
            results.append((self.arax_drug_curie_list[drug_idx][0],self.arax_drug_curie_list[drug_idx][1],np_list[drug_idx]))
            if(len(results)==n):break
        #results = [(arax_drug_curie_list[x][0],arax_drug_curie_list[x][1],np_list[x]) for x in ind]
        return results
    
    def yieldDrugsForMONDO(self,mondo,check_for_tps=False):
        drug_map = ARAX_Drug_Map(mondo,self)
        disease_idx = self.arax_mondo_to_idx[mondo]
        np_list = drug_map.arax_npz[drug_map.mondo_bin][drug_map.mondo_sub_idx]
        #https://stackoverflow.com/questions/6910641/how-do-i-get-indices-of-n-maximum-values-in-a-numpy-array
        for drug_idx in range(np_list.size):
            drug_curie = self.arax_drug_curie_list[drug_idx][0]
            drug_name =  self.arax_drug_curie_list[drug_idx][1]
            drug_mondo_score = np_list[drug_idx]
            if(check_for_tps):
                if((disease_idx,drug_idx) in tp_pairs):continue
            yield (drug_curie,drug_name,drug_mondo_score)
            
        #results = [(arax_drug_curie_list[x][0],arax_drug_curie_list[x][1],np_list[x]) for x in ind]
        return



#This class builds the second layer of the dictonary. So if x=ARAX_DICT["MONDO:1234567"],
# then when you call x["PUBCHEM:8910"], we are using this code.
# Basically this lets us map DrugNames to their appropirate identifiers in the numpy arrays.
class ARAX_Drug_Map(Mapping):
    def __init__(self,  mondo_idx, mondo_map):
        #If the provided MONDO is an integer, we assume someone already converted
        #the MONDO to the index in the big numpy matrix. Otherwise we assume mondo_idx 
        # in the format "MONDO:1234567" and look up it's index in our big dictonary.
        self.mondo_map = mondo_map
        if(type(mondo_idx)!=int):
            if(mondo_idx in self.mondo_map.arax_mondo_to_idx):
                mondo_idx = self.mondo_map.arax_mondo_to_idx[mondo_idx]
            else: raise ValueError
        self.mondo_idx = mondo_idx
        self.mondo_bin = mondo_idx // 1000
        self.mondo_sub_idx = mondo_idx % 1000
        self.arax_npz = self.mondo_map.arax_npz
    def __getitem__(self, key):
        if(type(key)!=int):
            if(key in self.mondo_map.arax_drug_curie_to_idx):
                key = self.mondo_map.arax_drug_curie_to_idx[key]
            else: raise ValueError(f"{key} cannot be found in our collection of valid drug identifiers.")
        return self.arax_npz[self.mondo_bin][self.mondo_sub_idx][key]
    def __iter__(self):
        return iter(self.mondo_map.arax_drug_curie_list)
    def __len__(self):
        return len(self.mondo_map.arax_drug_curie_to_idx)

def getIterForEachDrugForMONDO(mondo):
    yield "x"

def getScoreAndPercentile(mondo,chem):
#    CHEMBL.COMPOUND:CHEMBL1743073
    arax_dict = ARAX_MONDO_Map()
#    score = arax_dict[mondo][chem]
    a = arax_dict.getTopNHitsForMONDOFromArax(mondo,100000,False)
    #print(len(a))
    for i,h in enumerate(a):
        if(h[0]==chem):return f"{mondo}/{h[0]} Score:{h[2]:.4f}; Percentile {i}/{len(a)}:{i/len(a):.4f}"

def getScoreAndPercentileTuples(mondo,chem_set,return_all=False):
    from get_mondo_for_ct import buildMONDOLabelsDict

#    CHEMBL.COMPOUND:CHEMBL1743073
    arax_dict = ARAX_MONDO_Map()
    mondo_label_dict = buildMONDOLabelsDict()
    mondo_label_dict = {mondo_label_dict[x] : x for x in mondo_label_dict}
#    score = arax_dict[mondo][chem]
    a = arax_dict.getTopNHitsForMONDOFromArax(mondo,100000,False)
    #print(len(a))
    
    for i,h in enumerate(a):
        if(return_all or h[0] in chem_set):
            chem_curie = h[0]
            chem_name = h[1]
            mondo_curie = mondo
            mondo_name = mondo_label_dict[mondo_curie]
            pair_score = h[2]
            percentile = (i+1)/len(a)
            rank = i+1
            yield((mondo_curie,mondo_name,chem_curie,chem_name,pair_score,rank,percentile))
            #return f"{mondo}/{h[0]} Score:{h[2]}; Percentile {i}/{len(a)}:{i/len(a):.4f}"

def runXDTDForSingleMONDOAgainstEveryDrug(mondo_curie):
    import joblib
    from get_mondo_for_ct import buildMONDOLabelsDict
    label_to_mondo = buildMONDOLabelsDict()
    mondo_label_dict = {label_to_mondo[x] : x for x in label_to_mondo}
    mondo_name = mondo_label_dict[mondo_curie]

    #Load all of the embeddings generated for input into the xDTD model.
    drug_info_file = "/home/ubuntu/ARAX/drug_info.txt"
    mondo_embs = np.load('/home/ubuntu/ARAX/KGML-xDTD/mondo_embeddings.npz')["arr_0"]
    drug_embeddings = np.load('/home/ubuntu/ARAX/KGML-xDTD/drug_embeddings.npz')["arr_0"]
    #Load the xDTD model.
    model_path = "/home/ubuntu/ARAX/KGML-xDTD/model_evaluation/models/kgml_xdtd/drp_module/model.pt"
    fitModel = joblib.load(model_path)


    #This block of code is to get the embedding for the specific MONDO identifier from xDTD
    (arax_mondo_to_idx,arax_mondo_list) = buildMondoList()
    arax_disease_num = arax_mondo_to_idx[mondo_curie]
    mondo_emb = mondo_embs[arax_disease_num]

    #This builds up a matrix of DISEASE_EMBEDDING/DRUG_EMBEDDING for input to the model.
    X1 = np.tile(mondo_emb, (drug_embeddings.shape[0], 1))
    X2 = drug_embeddings
    X = np.hstack((X2,X1))

    #The model predicts 3 columns for each drug. The middle column is the one we care about;
    #it's the drug-disease treatment prediction
    res_from_xdtd = fitModel.predict_proba(X)
    #print(res_temp.shape)
    results = []
    with open(drug_info_file) as f:
        next(f)#we need to skip the header
        for line in f:
            chem_curie, chem_name, index = line.split('\t')[1]
            pair_score = res_from_xdtd[index][1]
            results.append((pair_score, mondo_curie,mondo_name,chem_curie,chem_name))
    for i,(pair_score, mondo_curie,mondo_name,chem_curie,chem_name) in enumerate(sorted(results,reverse=True)):
        percentile = (i+1)/len(results)
        rank = i + 1
        yield((mondo_curie,mondo_name,chem_curie,chem_name,pair_score,rank,percentile))
        
#chordoma_all_results = ARAX_NPZ.runXDTDForSingleMONDOAgainstEveryDrug("MONDO:0008978")

def normDrugDisease(df):
    df = normalize_nodes.normalize_data_frame(df,"Disease")
    df = normalize_nodes.normalize_data_frame(df,"Drug")
    return df

def getARAXScoresDF(mondo_curie,all_possible_drugs=False,normalize_drugs=False):
    if(all_possible_drugs):
        score_iter = runXDTDForSingleMONDOAgainstEveryDrug(mondo_curie)
    else:
        score_iter = getScoreAndPercentileTuples(mondo_curie,set(),True)
    df = pd.DataFrame(score_iter,columns=["Disease_Idx",'Disease_Name','Drug_Idx',"Drug_Name",'Score',"Rank","Percentile"])
    if(normalize_drugs):
        df = normalize_nodes.normalize_data_frame(df,"Drug_Idx")
    return df


if(__name__=="__main__"):
    runXDTDForSingleMONDOAgainstEveryDrug("MONDO:0000004")
    df = getARAXScoresDF("MONDO:0000004")
    chordoma_all_results = runXDTDForSingleMONDOAgainstEveryDrug("MONDO:0008978")
    exit()
    arax_dict = ARAX_MONDO_Map()
    print(arax_dict["MONDO:0000004"]["MESH:D000316"])
    a = arax_dict.getTopNHitsForMONDOFromArax("MONDO:0015564",15)
    print(a,len(a))
    
#arax_mondo = ARAX_MONDO_Map()
#large_pairs = []
#for mondo in tqdm(arax_mondo):
#    ARAX_DICT is an instance
#    arax_dict = arax_mondo[mondo]
#    np_list = arax_dict.arax_npz[arax_dict.mondo_bin][arax_dict.mondo_sub_idx]
#    for i,chem in enumerate(arax_dict):
#        print(mondo,chem)
#        worst_top_1000_score = top_pairs[num_pairs-1][2]
        #score = arax_dict[mondo].arax_npz[ara]
#        score = np_list[i]
        
#        if(score > 0.97):
#            large_pairs.append((mondo,chem[0],score))
        #if(score>worst_top_1000_score):
        #    updateList(mondo,chem[0],score)
#       
#        print(mondo,chem,score)