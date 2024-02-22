import multiprocessing as mp
import csv
import time
import orjson
import os
import tqdm
import sys
#Lots of the multiprocessing in this code is based on this 
#SO answer https://stackoverflow.com/a/13530258

def read_log(mondo_json_file, q):
    #print(mondo_json_file)
    if("MONDO" not in mondo_json_file): return
    mondo_idx = mondo_json_file.split("MONDO:")[1].split(".json")[0]
#    mondo_idx = int(mondo_idx) 
    mondo_idx = "MONDO:" + mondo_idx
    with open(mondo_json_file) as f:
       f_str = f.read().replace("Infinity","null")
       j = orjson.loads(f_str)
    if(j['message'].get('results',None)==None):return

    for result in j['message']['results']:
        try:
           score = result['score']
        except:
           print(mondo_idx)
           return
        if(score==None):
            score=-1.0
        else: score = float(score)
        alt_mondo_idx = result['node_bindings']['disease'][0]['id']
    
        mondo_name = j['message']['knowledge_graph']['nodes'][alt_mondo_idx]['name']
        chemical_idx = result['node_bindings']['chemical'][0]['id']
        chemical_name = j['message']['knowledge_graph']['nodes'][chemical_idx]['name']
        q.put((mondo_idx,alt_mondo_idx,mondo_name,chemical_idx,chemical_name,score))
        #q.put((mondo_idx,chemical,score))
    return

#First argument becomes file name
if(len(sys.argv)>1): fn=sys.argv[1]
else:fn = "aragorn_results.csv"

def listener(q):
    '''listens for messages on the q, writes to file. '''

    with open(fn, 'w') as f:
        writer = csv.writer(f)
        writer.writerow(["mondo_idx","alt_mondo_idx","mondo_name","chemical_idx","chemical_name","score"])
        while 1:
            m = q.get()
            if m == 'kill':
                break
            #(mondo_idx,mondo_name,chemical_idx,chemical_name,score) = m
            (mondo_idx,alt_mondo_idx,mondo_name,chemical_idx,chemical_name,score) = m
            #f.write(f"{mondo_idx},{chemical_idx},{score}\n")
            writer.writerow([mondo_idx,alt_mondo_idx,mondo_name,chemical_idx,chemical_name,score])
            #f.write(f"{mondo_idx},{alt_mondo_idx},{mondo_name},{chemical_idx},{chemical_name},{score}\n")
            f.flush()



l = []

#RESULT_DIR="../aragorn_results/"
#cnt=0
    
#    (score,chem) = read_log("results/" + f)
#    if(chem!=None): 
#        chem = chem[0]["id"]
#    disease = f.strip(".json")
#    f2.write(str(score) + ',' + str(disease) + "," + str(chem))

def handleJobs():
    for job in tqdm.tqdm(jobs): 
        job.get()

if(__name__=="__main__"):
    manager = mp.Manager()
    q = manager.Queue()
    pool = mp.Pool(mp.cpu_count() + 2)

    watcher = pool.apply_async(listener, (q,))
   
    jobs = []
    with open("query_info.csv") as f:
        next(f) #header
        for line in f:
            l = line.strip().split(',')
            mondo_idx = l[0]
            results = int(l[4])
            if(results==-1):
                pass
                #print(mondo_idx)
            if(results>0):
                json_file_path=os.path.join("MONDO",mondo_idx+'.json')
                job = pool.apply_async(read_log, (json_file_path, q))
                jobs.append(job)
                break
                if(len(jobs)>1000):
                    handleJobs()
                    jobs = []
            
#    for json_file in tqdm.tqdm(os.listdir(RESULT_DIR)):
#        json_file_path=os.path.join(RESULT_DIR,json_file)
#        job = pool.apply_async(read_log, (json_file_path, q))
#        jobs.append(job)


    #now we are done, kill the listener
    time.sleep(5)
    q.put('kill')
    handleJobs()
    pool.close()
    pool.join()    
