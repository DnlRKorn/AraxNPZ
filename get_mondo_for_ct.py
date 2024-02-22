import clinical_trial_scripts as cts
import drug_central_sql
from collections.abc import Iterable
from collections import defaultdict
import csv



def buildMONDOExactSynDict(exclude_ambigious:bool = False) -> (dict,dict) :
    mondo_matches = defaultdict(set)
    direct_mondos = {}

    with open("/home/ubuntu/PYTHON_SCRIPTS/DATA/mondo_syn.csv") as f:
        reader = csv.reader(f)
        next(reader)
        for row in reader:
            mondo_url = row[0]
            mondo_idx = "MONDO:" + mondo_url.split("_")[1]
            term = ','.join(row[1:]).lower()
            mondo_matches[term].add(mondo_idx)
            direct_mondos[term] = mondo_idx
    if(exclude_ambigious):
        direct_mondos = {}
        ambigious_mondos = {}
        for term in mondo_matches:
            mondo_id_set = mondo_matches[term]
            if(len(mondo_id_set)==1):
                direct_mondos[term] = mondo_id_set.pop()
            else:
                ambigious_mondos[term] = mondo_id_set
    else:
        ambigious_mondos = {}
        for term in mondo_matches:
            mondo_id_set = mondo_matches[term]
            if(len(mondo_id_set)==1): 
                pass
            else:
                ambigious_mondos[term] = mondo_id_set
    return (direct_mondos,ambigious_mondos)

def buildMONDOSynDict(fname : str, exclude_ambigious:bool = False, debug:bool=False) -> dict:
    syn_to_mondo = {}

    with open(fname) as f:
        reader = csv.reader(f)
        next(reader)
        for row in reader:
            mondo_url = row[0]
            mondo_idx = "MONDO:" + mondo_url.split("_")[1]
            term = ','.join(row[1:]).lower()
            syn_to_mondo[term] = mondo_idx
    
    if(exclude_ambigious):
        ambigious = set()
        with open(fname) as f:
            reader = csv.reader(f)
            next(reader)
            for row in reader:
                mondo_url = row[0]
                mondo_idx = "MONDO:" + mondo_url.split("_")[1]
                term = ','.join(row[1:]).lower()
                hit = syn_to_mondo[term]
                if(hit!=mondo_idx):ambigious.add(term)
        for term in ambigious:
            syn_to_mondo.pop(term)
    return syn_to_mondo

def buildMONDOLabelsDict(exclude_ambigious:bool=False,debug:bool=False) -> dict:
    return buildMONDOSynDict("/home/ubuntu/PYTHON_SCRIPTS/DATA/mondo_labels.csv",exclude_ambigious,debug)

def buildMONDOExactSynDict(exclude_ambigious:bool=False,debug:bool=False) -> dict:
    return buildMONDOSynDict("/home/ubuntu/PYTHON_SCRIPTS/DATA/mondo_exact_synonyms.csv",exclude_ambigious,debug)

def buildMONDONarrowSynDict(exclude_ambigious:bool=False,debug:bool=False) -> dict:
    return buildMONDOSynDict("/home/ubuntu/PYTHON_SCRIPTS/DATA/mondo_narrow_synonyms.csv",exclude_ambigious,debug)

def buildMONDORelatedSynDict(exclude_ambigious:bool=False,debug:bool=False) -> dict:
    return buildMONDOSynDict("/home/ubuntu/PYTHON_SCRIPTS/DATA/mondo_related_synonyms.csv",exclude_ambigious,debug)

def buildMONDOBroadSynDict(exclude_ambigious:bool=False,debug:bool=False) -> dict:
    return buildMONDOSynDict("/home/ubuntu/PYTHON_SCRIPTS/DATA/mondo_broad_synonyms.csv",exclude_ambigious,debug)

def getMONDOIdxDict() -> (dict,dict):
    mondo_idx_to_name = buildMONDOSynDict("/home/ubuntu/PYTHON_SCRIPTS/DATA/mondo_labels.csv",True,False)
    mondo_name_to_idx = {}
    for (mondo_idx,mondo_name) in mondo_idx_to_name.items():
        mondo_name_to_idx[mondo_name] = mondo_idx
    return (mondo_name_to_idx,mondo_idx_to_name)

def getConditionList(nct_list: Iterable ) -> (set,dict) :
    condition_set = set()
    condition_to_mesh = {}
    condition_len_list = []
    for nct in nct_list:
        path = cts.getNCTFilePath(nct)
        root = cts.getRootFromPath(path)
        #root = root.xpath("//clinical_study/condition_browse/mesh_term")
        #return
        cond_list = cts.getPropertyList(root, "//clinical_study/condition",exclude_NA=True,getTextOnly=True)
        condition_len_list.append(len(cond_list))
        #print(nct,list(cond_list))
        mesh_list_for_trial = cts.getPropertyList(root,"//clinical_study/condition_browse/mesh_term",exclude_NA=True,getTextOnly=True)
        for cond_name in cond_list:
            cond_name = cond_name.lower()
            condition_set.add(cond_name)
            mesh_hit = None
            for mesh in mesh_list_for_trial:
                #Extremely simple search, exact MeSH term in indication name.
                #A lot of indications are things like "ibuprofen Formulation 1" which contains the word "ibuprofen"
                if(mesh.lower() in cond_name):
                    mesh_hit = mesh.lower()
                if(mesh_hit!=None):
                    condition_to_mesh[cond_name]=mesh_hit
    #return condition_len_list
    return (condition_set, condition_to_mesh)

#This step is very simple. All it does is look at the exact primary name
# given from MONDO and sees if we can get a hit off of that.
def findMONDOsFromDict(cond_iter : Iterable, matches_found : set, mondo_dict : dict) -> (dict,dict):
    condname_to_matchname = {}
    condname_to_mondo = {}
    for cond_name in cond_iter:
        #We check only those conditions for which we haven't yet found a match.
        if(cond_name not in matches_found):
            if(cond_name in mondo_dict):
                condname_to_matchname[cond_name] = cond_name
                condname_to_mondo[cond_name] = mondo_dict[cond_name]
    return (condname_to_matchname,condname_to_mondo)

def checkSubstringInMONDODict(cond_iter : Iterable, matches_found : set, mondo_dict : dict, step_num : int=-1) -> (dict,dict):
#    condname_to_matchname = {}
#    condname_to_mondo = {}
    #cond_name_to_inexact_matches = {}
    inexact_hits = []
    for cond_name in cond_iter:
        #We check only those conditions for which we haven't yet found a match.
        if(cond_name not in matches_found):
            cond_hits = [(cond_name,key,mondo_dict[key],step_num) for key in mondo_dict.keys() if key in cond_name and key!='disease']
            inexact_hits.extend(cond_hits)
            #cond_name_to_inexact_matches[cond_name] = hits
            #if(len(hit_names)>1):
                #print(cond_name,len(hit_names))
                #break
            '''if(len(hit_names)>=1):
                hit_name = hit_names[0]
                if(hit_name=='d'):
                    print(hit_names,cond_name,list(mondo_dict)[0:3])
                    print('d' in mondo_dict)
                    print(mondo_dict['d'])
                    exit()
                condname_to_matchname[cond_name] = hit_name
                condname_to_mondo[cond_name] = mondo_dict[hit_name]'''
    return inexact_hits
    #return cond_name_to_inexact_matches

def checkSubstringInCondNames(cond_iter : Iterable, matches_found : set, mondo_dict : dict, step_num : int=-1) -> (dict,dict):
#    condname_to_matchname = {}
#    condname_to_mondo = {}
#    condname_to_inexact_matches = {}
    inexact_hits = []
    for cond_name in cond_iter:
        #We check only those conditions for which we haven't yet found a match.
        if(cond_name not in matches_found):
            hits = [(cond_name,key,mondo_dict[key],step_num) for key in mondo_dict.keys() if cond_name in key]
            inexact_hits.extend(hits)
            #condname_to_inexact_matches[cond_name] = hits
            '''if(len(hit_names)>=1):
                hit_name = hit_names[0]
                condname_to_matchname[cond_name] = hit_name
                condname_to_mondo[cond_name] = mondo_dict[hit_name]'''
    return inexact_hits
    #return condname_to_inexact_matches

def print5(d,n=5):
    for i,key in enumerate(d):
        if(i==0):print("cond_name:match_name")
        if(i>=n):break
        print(f"{key}:{d[key]}")
    return



def getConditionMatches(condition_set : Iterable[str], debug : bool=False, check_condition_in_mondo : bool=True, check_mondo_in_condition : bool=True) -> (dict,dict,list):
    condname_to_matchname = {}
    condname_to_mondo = {}
    exact_matches_found = set()
    #condition_set, condition_to_mesh = getConditionList(nct_iter)
    if(debug):print(f"Step 0: We have {len(condition_set)} conditions.")
    debug_print_cnt = 3
    inexact_matches = []

    def updateOneToOneBasedOnDict(dict,exact_matches_found,step_num):
        tmp_condname_to_matchname,tmp_condname_to_mondo = findMONDOsFromDict(condition_set, exact_matches_found, dict)
        condname_to_matchname.update(tmp_condname_to_matchname)
        condname_to_mondo.update(tmp_condname_to_mondo)
        exact_matches_found.update(condname_to_mondo.keys())
        if(debug):print(f"Step {step_num}: We have found {len(condname_to_matchname)} hits. Ratio {len(condname_to_matchname)/len(condition_set)}")
        if(debug):print5(tmp_condname_to_matchname,debug_print_cnt)
        return exact_matches_found


    #
    def updateConditionContainedInMONDO(dict,exact_matches_found,step_num):
        tmp_inexact_matches = checkSubstringInMONDODict(condition_set, exact_matches_found, dict,step_num)
        inexact_matches.extend(tmp_inexact_matches)
        #condname_to_matchname.update(tmp_condname_to_matchname)
        #condname_to_mondo.update(tmp_condname_to_mondo)
        #exact_matches_found.update(condname_to_mondo.keys())
        if(debug):print(f"Step {step_num}: We have found {len(condname_to_matchname)} hits. Ratio {len(condname_to_matchname)/len(condition_set)}")
        if(debug):print5(tmp_condname_to_matchname,debug_print_cnt)
        return exact_matches_found

    def updateMONDOContainedInConditionName(dict,exact_matches_found,step_num):
        tmp_inexact_matches = checkSubstringInCondNames(condition_set, exact_matches_found, dict,step_num)
        inexact_matches.extend(tmp_inexact_matches)

        #condname_to_matchname.update(tmp_condname_to_matchname)
        #condname_to_mondo.update(tmp_condname_to_mondo)
        #exact_matches_found.update(condname_to_mondo.keys())
        if(debug):print(f"Step {step_num}: We have found {len(condname_to_matchname)} hits. Ratio {len(condname_to_matchname)/len(condition_set)}")
        if(debug):print5(tmp_condname_to_matchname,debug_print_cnt)
        return exact_matches_found
    #Step 1: Direct matches from MONDO names.
    #Turns out the MONDO name list is corrupted; will need to fix later.
    #mondo_labels = buildMONDOLabelsDict(True) 
    mondo_idx_to_name,mondo_name_to_idx = getMONDOIdxDict()
    mondo_labels = buildMONDOLabelsDict(True)
    exact_matches_found = updateOneToOneBasedOnDict(mondo_labels,exact_matches_found,1)

    #Step 2: Direct matches from unambigious synonyms.
    exact_syns = buildMONDOExactSynDict(True)
    exact_matches_found = updateOneToOneBasedOnDict(exact_syns,exact_matches_found,2)

    #Step 3
    narrow_syns = buildMONDONarrowSynDict(True)
    exact_matches_found = updateOneToOneBasedOnDict(narrow_syns,exact_matches_found,3)

    #Step 4
    related_syns = buildMONDORelatedSynDict(True)
    exact_matches_found = updateOneToOneBasedOnDict(related_syns,exact_matches_found,4)
    
    #Step 5
    broad_syns = buildMONDOBroadSynDict(True)
    exact_matches_found = updateOneToOneBasedOnDict(broad_syns,exact_matches_found,5)

    if(check_condition_in_mondo):
        matches_found = updateConditionContainedInMONDO(mondo_labels,exact_matches_found,6)
        matches_found = updateConditionContainedInMONDO(exact_syns,exact_matches_found,7)
        matches_found = updateConditionContainedInMONDO(narrow_syns,exact_matches_found,8)
        matches_found = updateConditionContainedInMONDO(related_syns,exact_matches_found,9)
        matches_found = updateConditionContainedInMONDO(broad_syns,exact_matches_found,10)

    if(check_mondo_in_condition):
        matches_found = updateMONDOContainedInConditionName(mondo_labels,exact_matches_found,11)
        matches_found = updateMONDOContainedInConditionName(exact_syns,exact_matches_found,12)
        matches_found = updateMONDOContainedInConditionName(narrow_syns,exact_matches_found,13)
        matches_found = updateMONDOContainedInConditionName(related_syns,exact_matches_found,14)
        matches_found = updateMONDOContainedInConditionName(broad_syns,exact_matches_found,15)

    #inexact matches is as follows(cond_name,key,mondo_dict[key],step_num) 
    return (condname_to_matchname,condname_to_mondo,inexact_matches)

def getHits(nct_iter : Iterable[str], debug : bool=False) -> (dict,dict):
    condition_set, condition_to_mesh = getConditionList(nct_iter)
    (condname_to_matchname,condname_to_mondo,inexact_matches) = getConditionMatches(condition_set, debug)
    return (condname_to_matchname,condname_to_mondo,inexact_matches)

def getOneCondition(condition: str, debug:bool=False):
    condition = condition.lower()
    (condname_to_matchname,condname_to_mondo,inexact_matches) = getConditionMatches([condition], debug)
    return (condname_to_matchname.get(condition,None),condname_to_mondo.get(condition,None),inexact_matches)

def trialGen():
    with open("/home/ubuntu/PYTHON_SCRIPTS/DATA/first_post_after_march.txt") as f:
        for line in f:
            nct = line.strip()
            yield nct

if(__name__=="__main__"):
    print(getConditionMatches(["rhabdomyolysis"],True))
    
#    mondo_labels = buildMONDOLabelsDict(True)
#    print(checkSubstringInMONDODict(["breast cancer metastatic"], set(), mondo_labels))
#    print(checkSubstringInCondNames(["breast cancer metastatic"], set(), mondo_labels))
    #print(tmp_condname_to_matchname)
#    buildMONDORelatedSynDict()
#    nct_list = []
#    with open("/home/ubuntu/PYTHON_SCRIPTS/DATA/p-value-cnt.txt") as f:
#        next(f) #header
#        for line in f:
#            nct = line.split('.')[0]

            #(fname, pval) = line.split(":")
            #path = cts.getNCTFilePath(fname)
            #root = cts.getRootFromPath(path)

            #if(not cts.check005LessPVal(root)):continue
#            nct_list.append(nct)
#            if(len(nct_list)>=100):break
    #print("NCT_LIST",len(nct_list))
#    len_list = getConditionList(nct_list)
    #getMONDOIdxDict()
    print(len(list(trialGen())))
#    print(len(getConditionList(trialGen())[0]))
    #len_list = getConditionList(trialGen())

    #import pandas as pd
    #sr = pd.Series(len_list)
    #df = sr.to_frame()
    #df2 = df.value_counts()
    #df2.to_csv("/home/ubuntu/PYTHON_SCRIPTS/DATA/ct_condition_counts.csv")
#    print(df2)
    
    #(condname_to_matchname,condname_to_mondo,inexact_matches) = getHits(list(trialGen())[0:20],False)
    #print(condname_to_matchname)
    #print(inexact_matches)
