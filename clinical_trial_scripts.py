#import xml.etree.ElementTree as ET
import itertools
import lxml.etree as ET
import os
from typing import List
import csv
#import dateutil]
from datetime import datetime
from dateutil.parser import parse
from collections.abc import Iterable
#from tqdm.notebook import tqdm
from collections import defaultdict

#import drug_central_sql
#import importlib
#importlib.reload(drug_central_sql)

#https://www.w3schools.com/xml/xpath_syntax.asp

#tree = ET.parse('country_data.xml')
#root = tree.getroot()

#In the AllPublicXML release from clinicaltrials.gov, each clinical trial is
# stored in a directory NCT1234xxxx. So to get NCT12345678.xml we need the
# file path NCT1234xxxx/NCT12345678.xml . All this code does is construct
# that directory header and append it to the front of the file path.
def getNCTFilePath(nct: str):
    if(not nct.endswith(".xml")):nct = nct+".xml"
    clinical_trial_dir = "/home/ubuntu/PYTHON_SCRIPTS/DATA/AllPublicXML"
    nct_directory = nct[0:7] + "xxxx"
    path = os.path.join(clinical_trial_dir,nct_directory,nct)
    return path

#Opens the XML file passed in through the "path" parameter. 
def getRootFromPath(path: str) -> ET._Element:
    #First see if the path is valid, if it isn't, try to complete the path.
    #If that doesn't work, return None.
    if(not os.path.isfile(path)): 
        path = getNCTFilePath(path)
        if(not os.path.isfile(path)): return None
    tree = ET.parse(path)
    root = tree.getroot()
    return root

#Finds the first matching item with the tag "prop" and
# returns the text from that element.
def getPropertyFromRoot(root: ET._Element, prop: str) -> str:
    x = root.find(prop)
    if(x==None):return "N/A"
    return x.text

def getElementFromXPath(root: ET._Element, xpath: str) -> List[ET._Element]:
    x = root.xpath(xpath)
    return x

#Searches
# the XML for a "condition" element immediately off of the first level
# of the XML tree. If it cannot find it, we return N/A. If we find it
# return the text (the condition).
def getConditionFromRoot(root: ET._Element) -> str:
    y = getPropertyFromRoot(root, "condition")
    return y

def getPhaseFromRoot(root) -> str:
    y = getPropertyFromRoot(root, "phase")
    return y

def getLink(NCT) -> str:
    link = "https://clinicaltrials.gov/ct2/show/"
    return link + NCT

def getPVal_iter(root) -> Iterable:
    return root.iter("p_value")

#If we can find a single p_value field with value <0.05, return True immediately.
def check005LessPVal(root) -> bool:
    pval_iter = getPVal_iter(root)
    for x in find005LessPValIter(pval_iter):
        return True
    return False

#Builds a python iterator object. Goes through the p_value's from the XML
# tree one at a time. If the value of the text is <=0.05, yield it.
def find005LessPValIter(p_val_iter) -> Iterable:
    for p_val in p_val_iter:
        outcome_type = p_val.xpath("../../../type")[0].text
        
        if(outcome_type!="Primary"):continue
        
        ptext = p_val.text

        ptext = ptext.replace(" ","")
        ptext = ptext.replace("=","")
        ptext = ptext.replace("p","")
        ptext = ptext.replace(",",".")
#       print(trial_name,ptext)

        if(ptext[0]==">"):
            p_float = float(ptext.replace(">",""))
            #We found a good probablility, stop searching
            if(p_float < 0.05):
                yield p_val
        elif("<" in ptext):
            p_float = float(ptext.replace("<",""))
            #We found a good probablility, stop searching
            if(p_float <= 0.05):
                yield p_val
        else:
            p_float = float(ptext)
            #We found a good probablility, stop searching.
            if(p_float <= 0.05):
                yield p_val
    return False

#Clinical trials can have a lot of weirdness in how p-values are written
# that make pythons float interperater unhappy. The main things that cause us
# issues are random spaces and putting the letter "p" at the beginning of 
# p-values.
def pvalueToFloat(ptext):
    ptext = ptext.replace(" ","")
    ptext = ptext.replace("=","")
    ptext = ptext.replace("p","")
    ptext = ptext.replace(",",".")
    ptext = ptext.replace(">","")
    ptext = ptext.replace("<","")
    return float(ptext)
    
def getEnrollmentAndType(root):
    x = root.find("enrollment")
    if(x==None): return ("N/A","N/A")
    typ = x.get("type")
    text = x.text
    return (typ,text)

def getIntervention(root):
    inds = root.xpath("//clinical_study/intervention")
    drug_cnt = 0
    bio_cnt = 0
    ind_info = []
    ind_list = []
    for ind in inds:
        ind_type = getPropertyString(ind,"intervention_type")
        ind_name = getPropertyString(ind,"intervention_name")
        ind_description = getPropertyString(ind,"description") 
        if(ind_type=="Drug"): drug_cnt+=1
        if(ind_type=="Biological"): bio_cnt+=1
        ind_info.append(f"{ind_type}:{ind_name}/{ind_description}")
        ind_list.append((ind_type,ind_name,ind_description))
    int_info_str = "|".join(ind_info)
    return (int_info_str,drug_cnt,bio_cnt, ind_list)
        

#Iterates through all XML files we have saved from clincal trials
# and yields an LXML object which parses that XML.
def getTrialIter() -> Iterable[ET._Element]:
    with open("p-value-cnt.txt") as f:
        next(f) #header
        for line in f:
            (fname, _) = line.split(":")
            path = getNCTFilePath(fname)
            root = getRootFromPath(path)
            yield root
    
def getTrialIter005LessPVal() -> Iterable[ET._Element]:
    for root in getTrialIter():
        if(check005LessPVal(root)): yield root
            
def getGroupInfo(outcome):
    group_info = []
    groups = outcome.xpath("measure/analyzed_list/analyzed/count_list/count")
    for group in groups:
        group_id = group.get("group_id")
        value = group.get("value")
        group_info.append(group_id + "/" + value)
    return "|".join(group_info)

def getPropertyString(node,property_xpath,exclude_NA=False) -> str:
    vals = node.xpath(property_xpath)
    #XPath returns a list of all xml nodes matching our pattern. What
    # we do is grab the text from each node, and combine them in a list.
    if(exclude_NA): vals_text_list = [val.text for val in vals if val.text!="NA"]
    else:vals_text_list = [val.text for val in vals]
            #Combines all text fields into a single value. So
            # 3 fields would be stored in a list [0.01,0.02,0.04]
            # becomes 0.01|0.02|0.04
    if(len(vals_text_list)==0):return "N/A"
    property_string = "|".join(vals_text_list)
    return property_string

def getPropertyList(node : ET._Element,property_xpath : str,exclude_NA : bool =False,getTextOnly : bool =False) -> List:
    vals = node.xpath(property_xpath)
    #XPath returns a list of all xml nodes matching our pattern. What
    # we do is grab the text from each node, and combine them in a list.
    if(exclude_NA): val_list = [val for val in vals if val.text!="NA"]
    else: val_list = [val for val in vals]
    
    if(len(val_list)==0):return []
    
    if(getTextOnly): return [val.text for val in val_list]
    else: return val_list

def iterClinicalTrials() -> Iterable[str]:
    clinical_trial_dir = "/home/ubuntu/PYTHON_SCRIPTS/DATA/AllPublicXML"
    for sub_dir in os.listdir(clinical_trial_dir):
        if("NCT" not in sub_dir):continue
#        print(sub_dir)
        for trial_name in os.listdir(os.path.join(clinical_trial_dir,sub_dir)):
            yield trial_name
#len(list(iterClinicalTrials()))

#This adds in a check if the date is "Actual" or "Estimated"
def getDateProperty(root,prop,only_actual=True) -> str:
    x = root.find(prop)
    #if(x==None):return "N/A"
    
    #This will return true if the "Type" of the date
    # property is either non-existent or "Actual". This
    # is just how the ClinicalTrials xml format. It
    # may also be "Estimated" or "Anticipated". Which
    # we don't want.
    if(x!=None):
        if(only_actual):
            if(x.attrib.get("type","Actual")=="Actual"):
                return x.text
        else:
            return x.text
    return "N/A"

def getDatesFromNCT(nct : str, only_actual :bool =True):
    trial_path = getNCTFilePath(nct)
    #Opens the XML file passed in through the "path" parameter. 
    root = getRootFromPath(trial_path)

    start_date = getDateProperty(root,"start_date",only_actual)
    study_first_submitted = getDateProperty(root,"study_first_submitted",only_actual)
    completion_date = getDateProperty(root,"completion_date",only_actual)
    last_update_submitted = getDateProperty(root,"last_update_submitted",only_actual)
    results_first_posted = getDateProperty(root,"results_first_posted",only_actual)
    last_update_posted = getDateProperty(root,"last_update_posted",only_actual)
    verification_date = getDateProperty(root,"verification_date",only_actual)
    primary_completion_date = getDateProperty(root,"primary_completion_date",only_actual)
    phase = getPropertyFromRoot(root,"phase")
    #["start_date","study_first_submitted","completion_date","last_update_submitted","last_update_posted","verification_date"]
    d = {"start_date":start_date,"study_first_submitted":study_first_submitted,
         "completion_date":completion_date,"last_update_submitted":last_update_submitted,
         "last_update_posted":last_update_posted,"verification_date":verification_date,
         "primary_completion_date":primary_completion_date,"results_first_posted":results_first_posted}
#    return (start_date,study_first_submitted,completion_date,last_update_submitted,last_update_posted,verification_date)
    return d

def checkIfTrialIsPastDate(nct : str, target_date : str, key : str) -> bool:
#    f
    trial_path = getNCTFilePath(nct)
    root = getRootFromPath(trial_path)    
    date_str = getDateProperty(root,key)
#    date_dict = getDatesFromNCT(nct)
#    date_str = date_dict[key]
    if(date_str=="N/A"):return False
   
    date_obj = parse(date_str)
    if(date_obj > target_date): return True
    return False
#    for date_str in date_dict:
#        if(date_str=="N/A"):continue
##        date_obj = parse(date_str)

#This function is meant to help with multithreading.
# Basically it will only make any meaningful difference to 
# code speed if the list is from multiprocessing
def checkTrialFunc(nct : str,target_date : str,key : str, mp_list: List):
    if(checkIfTrialIsPastDate(nct,target_date,key)):
        #print(nct)
        mp_list.append(nct)
    return

def multiThreadCheckTrials(target_date,key,chunksize=1024,pbar=None):
    import multiprocessing as mp
    manager = mp.Manager()
    #q = manager.Queue()
    mp_list = manager.list()
    #pool = mp.Pool(mp.cpu_count() + 2)
 
    #watcher = pool.apply_async(listener, (q,))
    def funcArgumentGen(target_date,key,mp_list,pbar=None):
        #print('hi')
        #yield(("NCT05751629",target_date,key,l))
        #print('hi2')
        for _,nct in enumerate(iterClinicalTrials()):
        #for i, nct in [(0,"NCT03718637.xml")]:
            yield ((nct,target_date,key,mp_list))
            if(type(pbar)!=type(None)): pbar.update(1)
            #if(i>50000):return
    #print('hi')
    pbar.total = 4000
    pbar.update(4)
    trial_len = len(list(funcArgumentGen('a','b','c',None)))
    with mp.Pool(mp.cpu_count() + 2) as pool:
        if(type(pbar)!=type(None)):pbar.total=trial_len
        _ = pool.starmap_async(checkTrialFunc,funcArgumentGen(target_date,key,mp_list,pbar), chunksize=chunksize)
        pool.close()
        pool.join()
    return list(mp_list)
    #while(not q.empty()):
    #    l.append(q.get)
    #return l
import time

def buildTableForNCTs(nct_iter: Iterable, fname: str):
    from tqdm import tqdm
    fields = ['nct','url','trial_name','phase','conditions','interventions','start_date', 'study_first_submitted', 'completion_date', 'last_update_submitted', 'last_update_posted', 'verification_date', 'primary_completion_date', 'results_first_posted']
    with open(fname,'w') as f:
        d_writer = csv.DictWriter(f,fields)
        d_writer.writeheader()
        for nct in tqdm(nct_iter,total=500000):
            d = getDatesFromNCT(nct)
#            for k in d:
#                if(d[k]!="N/A"):
#                    d[k]=parse(d[k])
            trial_path = getNCTFilePath(nct)
            root = getRootFromPath(trial_path)
            nct = nct.replace(".xml",'')
            d['nct'] = nct
            d['url'] = "https://classic.clinicaltrials.gov/ct2/show/" + nct
            d['trial_name'] = getPropertyFromRoot(root,'brief_title')
            d['phase'] = getPropertyFromRoot(root,"phase")
            (int_info_str,drug_cnt,bio_cnt, ind_list) = getIntervention(root)
            d['interventions'] = int_info_str
            conditions = getPropertyString(root,"//clinical_study/condition",True)
            d['conditions'] = conditions
            d_writer.writerow(d)


    #return None

if(__name__=="__main__"):
    #<completion_date type="Anticipated">April 14, 2025</completion_date>
    #<primary_completion_date type="Anticipated">April 14, 2024</primary_completion_date>
    #print(getNCTFilePath("NCT05751629"))
    #exit()

    nct = "NCT00501761"
    d = getDatesFromNCT(nct)
    print(d.keys())
    buildTableForNCTs(iterClinicalTrials(),'/home/ubuntu/PYTHON_SCRIPTS/DATA/test.csv')
    exit()
    """     with open("/home/ubuntu/PYTHON_SCRIPTS/DATA/after_march_ct.csv") as f:
        for i,line in enumerate(f):
            nct = line.strip()
            trial_path = getNCTFilePath(nct)
            root = getRootFromPath(trial_path)    
            l = getElementFromXPath(root,"//clinical_study/completion_date")

#    date_str = getPropertyFromRoot(root,"completion_date")
#<completion_date type="Anticipated">April 14, 2025</completion_date>
#  <primary_completion_date type="Anticipated">April 14, 2024</primary_completion_date>
            b1 = l[0].attrib.get("type") not in ["Estimated","Anticipated"]
            b2 = l[0].attrib.get("type","Actual")=="Actual"
            if(b1!=b2):
                print(i,nct,l[0].attrib.get("type") not in ["Estimated","Anticipated"],l[0].attrib.get("type","Actual")=="Actual")
    #        if(i>200):break
    exit() 
    for i in range(1,8):
        chunksize = 4**i
        t1 = datetime.now()
        
        t2 = datetime.now()
        print(chunksize, (t2-t1).total_seconds())
    from tqdm import tqdm
    pbar = tqdm(total=55)
    chunksize=1024

    t1 = datetime.now()
    l = multiThreadCheckTrials("Feb 28, 2023","results_first_posted",chunksize,pbar)    
    t2 = datetime.now()
    print(chunksize,'async', (t2-t1).total_seconds())

    pbar.close()
    with open("/home/ubuntu/PYTHON_SCRIPTS/DATA/first_post_after_march.txt",'w') as f:
        for x in l:
            f.write(x + '\n')
"""
#    x = checkIfTrialIsPastDate("NCT05751629","Feb 28, 2023")
    
    #print(l)
    '''
    Trial time runs
    1 130.779406
10 30.232904

100 13.917001
1000 12.627759
10000 11.986119
100000 36.336085
1000000 39.613132
10000000 37.356472
100000000 40.133161
1000000000 39.924551
(base) ubuntu@ip-172-31-90-56:~$ /home/ubuntu/miniconda3/bin/python /home/ubuntu/PYTHON_SCRIPTS/clinical_trial_scripts.py
10 291.482965
100 273.987506
1000 298.173403
10000 323.129458'''