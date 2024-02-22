#!/usr/bin/env python
# coding: utf-8

# In[1]:


import psycopg2
from psycopg2.extras import execute_values


# In[ ]:





# In[2]:


connection = psycopg2.connect(user="drugman",
                              password="dosage",
                              host="unmtid-dbs.net",
                              port="5433",
                              database="drugcentral")

def testConn():
    # Create a cursor to perform database operations
    cursor = connection.cursor()
    # Print PostgreSQL details
    print("PostgreSQL server information")
    print(connection.get_dsn_parameters(), "\n")
    # Executing a SQL query
    cursor.execute("SELECT version();")
    # Fetch result
    record = cursor.fetchone()
    print("You are connected to - ", record, "\n")


# In[3]:


def getCursor():
    global connection
    try:
        if(connection.isolation_level==None or connection.closed!=0):
            connection = psycopg2.connect(user="drugman",
                              password="dosage",
                              host="unmtid-dbs.net",
                              port="5433",
                              database="drugcentral")    
    except OperationalError as oe:
        connection = psycopg2.connect(user="drugman",
                              password="dosage",
                              host="unmtid-dbs.net",
                              port="5433",
                              database="drugcentral")
    cursor = connection.cursor()
    return cursor


# ### Get off-patent drugs

# In[4]:


def offPatent():
    query = cursor.mogrify('''SELECT s.id, s.name
    FROM public.structures s 
    WHERE s.status='OFP'
    LIMIT 10;''')
    print(query)
    cursor = connection.cursor()

    cursor.execute(query)
    record = cursor.fetchone()
    print(record)


# ### Get drugs with Off-Label Uses

# In[5]:


def offLabel():
    query = cursor.mogrify('''SELECT DISTINCT s.id, s.name 
    FROM public.omop_relationship or2
    JOIN public.structures s on or2.struct_id=s.id
    WHERE or2.relationship_name='off-label use'
    AND or2.concept_name<>'Coronavirus infection'
    LIMIT 10;''')

    print(query)
    cursor = connection.cursor()

    cursor.execute(query)
    for record in cursor:
        print(record)


# In[6]:


def countOffLabel():
    #Count the number of offLabel drugs in DrugCentral.
    query = cursor.mogrify('''SELECT COUNT(DISTINCT(s.name)) 
    FROM public.omop_relationship or2
    JOIN public.structures s on or2.struct_id=s.id
    WHERE or2.relationship_name='off-label use'
    AND or2.concept_name<>'Coronavirus infection'
    LIMIT 10;''')

    print(query)
    cursor = connection.cursor()

    cursor.execute(query)
    for record in cursor:
        print(record)


# ### Get drugs and their off-label uses

# In[7]:


#This function gets drugs and the off label indications they serve.
def getDrugsWithOffLabel():
    query = cursor.mogrify('''SELECT s.id, s.name, or2.concept_id, or2.concept_name 
    FROM public.omop_relationship or2
    JOIN public.structures s on or2.struct_id=s.id
    WHERE or2.relationship_name='off-label use'
    AND or2.concept_name<>'Coronavirus infection'
    LIMIT 10;''')

    print(query)
    cursor = connection.cursor()

    cursor.execute(query)
    for record in cursor:
        print(record)


# ### Query with list parameter example

# In[8]:


def listParamTest():
    query = cursor.mogrify('''SELECT s.lname, or2.concept_name, or2.relationship_name 
     from public.synonyms s 
     join public.omop_relationship or2 on or2.struct_id=s.id
     where s.lname in %s;''', (("choline fenofibrate","test2","test3"),))#(("choline fenofibrate","TEST"),))
    print(query)
    cursor = connection.cursor()

    cursor.execute(query)
    record = cursor.fetchone()
    print(record)

    #cursor.execute("SELECT s.lname, or2.concept_name, or2.relationship_name  from public.synonyms s  join public.omop_relationship or2 on or2.struct_id=s.id where s.lname in ('choline fenofibrate');")
    #record = cursor.fetchone()
    #print(record)


# In[ ]:





# In[ ]:





# In[9]:


def checkForMeshId(mesh_list):
    cursor = getCursor()
    query = '''SELECT identifier, s."name", s.id 
    FROM public.identifier i
    join public.structures s on s.id=struct_id 
    where (id_type='MESH_DESCRIPTOR_UI' or id_type='MESH_SUPPLEMENTAL_RECORD_UI')
    and i.identifier in %s;'''
    query = cursor.mogrify(query,(mesh_list,))
    cursor = connection.cursor()
    cursor.execute(query)
    return list(cursor)


# In[10]:


mesh_ids = '''D004953
D002185
D000804
D000080045
D052203'''
mesh_ids = mesh_ids.split()
x = checkForMeshId(tuple(mesh_ids))


# In[11]:


#for a in x:
#    print(a)


# In[12]:


def checkSynonym(syn):
    cursor = getCursor()
    query = '''SELECT s.id, syn.lname, s."name"
    FROM public.synonyms syn
    JOIN public.structures s on s.id=syn.id
    WHERE syn.lname=%s;'''
    query = cursor.mogrify(query,(syn,))
    cursor = connection.cursor()
    cursor.execute(query)
    return list(cursor)


# In[13]:


def checkSynonymLIKE(syn):
    cursor = getCursor()
    query = '''SELECT DISTINCT(s.id)
    FROM public.synonyms syn
    JOIN public.structures s on s.id=syn.id
    WHERE syn.lname LIKE %s;'''
    query = cursor.mogrify(query,('%' + syn + '%',))
    #print(query)
    cursor = connection.cursor()
    cursor.execute(query)
    return list(cursor)


# In[14]:


#checkSynonym("sirolimus")


# In[15]:


#checkSynonymLIKE("beclomethasone")


# In[16]:


#https://stackoverflow.com/questions/65412161/execute-a-query-for-multiple-sets-of-parameters-with-psycopg2
'''
data = (
    (1, '2020-11-19'),
    (1, '2020-11-20'),
    (1, '2020-11-21'),
    (2, '2020-11-19'),
    (2, '2020-11-20'),
    (2, '2020-11-21')
)
        
query = """
    with data(song_id, date) as (
        values %s
    )
    select t.*
    from my_table t
    join data d 
    on t.song_id = d.song_id and t.date = d.date::date
"""
execute_values(cursor, query, data)
results = cursor.fetchall()
'''


# In[17]:


def bulkSynQuery(syn_list):
    query = """
    with data(lname) as (
        values %s
    )
    select syn.lname, syn.id
    from public.synonyms syn
    join data d 
    on syn.lname = d.lname
    """
    cursor = getCursor()
    #data = (("acyclovir",),("doxorubicin",),("albuterol",))
    syn_tuple_list = tuple([(x,) for x in syn_list])
    execute_values(cursor, query, syn_tuple_list)
    #results = cursor.fetchall()
    #return results
    return cursor

def bulkSynQueryIter(drug_list, range_val=50):
    query = """
    with data(lname) as (
        values %s
    )
    select syn.lname, syn.id
    from public.synonyms syn
    join data d 
    on syn.lname = d.lname
    """
    cursor = getCursor()
    #data = (("acyclovir",),("doxorubicin",),("albuterol",))
    #for i in tqdm(range(0, len(drug_list), range_val)):
    for i in range(0, len(drug_list), range_val):
        chunk = drug_list[i:i + range_val]
        chunk_tuple_list = tuple([(x,) for x in chunk])
        execute_values(cursor, query, chunk_tuple_list)
        for x in cursor:
            yield x


# In[18]:


drug_list = ["acyclovir","doxorubicin","albuterol"]


# In[19]:


#res = bulkSynQuery(drug_list)
#for x in res:print(x)


# In[20]:


#results = cursor.fetchall()


# In[21]:


#results


# In[22]:


def testBulk():
    cnt = 0
    for (drug_name,cd_id) in bulkSynQueryIter(drug_list * 100,100):
        print(drug_name)
        cnt+=1


# In[23]:


def bulkSynQuerySubstringIter(drug_list, range_val=50):
    query = """WITH data(lname) as (
    values %s 
    ) 
    SELECT d.lname, syn.id, syn.lname 
    FROM data d, public.synonyms syn 
    WHERE syn.lname LIKE CONCAT('%%',d.lname,'%%') """
    cursor = getCursor()
    #data = (("acyclovir",),("doxorubicin",),("albuterol",))
    #for i in tqdm(range(0, len(drug_list), range_val)):
    for i in range(0, len(drug_list), range_val):
        chunk = drug_list[i:i + range_val]
        chunk_tuple_list = tuple([(x,) for x in chunk])
        execute_values(cursor, query, chunk_tuple_list)
        for x in cursor:
            yield x
            
def bulkSynQuerySubstringReverseIter(drug_list, range_val=50):
    query = """WITH data(lname) as (
    values %s 
    ) 
    SELECT d.lname, syn.id, syn.lname 
    FROM data d, public.synonyms syn 
    WHERE d.lname LIKE CONCAT('%%',syn.lname,'%%') """
    cursor = getCursor()
    #data = (("acyclovir",),("doxorubicin",),("albuterol",))
    #for i in tqdm(range(0, len(drug_list), range_val)):
    for i in range(0, len(drug_list), range_val):
        chunk = drug_list[i:i + range_val]
        chunk_tuple_list = tuple([(x,) for x in chunk])
        execute_values(cursor, query, chunk_tuple_list)
        for x in cursor:
            yield x


# In[24]:


def testSubstring():
    cnt = 0
    drug_list = ["acyclovir","doxorubicin","albuterol"]
    drug_list = ["acy"]
    #for (drug_query,drug_name,cd_id) in bulkSynQuerySubstringIter(drug_list * 50,10):
    for x in bulkSynQuerySubstringIter(drug_list * 50,50):
        print(x)
        cnt+=1
    print(cnt)


# In[57]:


def checkForDBID(dbid):
    cursor = getCursor()
    query = '''SELECT identifier, struct_id
    FROM public.identifier
    WHERE id_type='DRUGBANK_ID'
    AND identifier=%s;'''
    query = cursor.mogrify(query,(dbid,))
    #print(query)
    cursor = connection.cursor()
    cursor.execute(query)
    return list(cursor)

def checkForDBIDBulk(dbid_list,range_val=50):
    cursor = getCursor()
    query = '''
    WITH data(dbid) as (
    values %s
    ) 
    SELECT identifier, struct_id
    FROM data d
    JOIN public.identifier p ON d.dbid=p.identifier
    WHERE p.id_type='DRUGBANK_ID';
    '''
    #query = cursor.mogrify(query,(dbid,))
    cursor = getCursor()
    for i in range(0, len(dbid_list), range_val):
        chunk = dbid_list[i:i + range_val]
        chunk_tuple_list = tuple([(x,) for x in chunk])
        execute_values(cursor, query, chunk_tuple_list)
        for x in cursor:
            yield x


# In[59]:


#list(checkForDBIDBulk(['DB00158','DB00514']))


# In[27]:


def bulkDCIDQueryIter(dcid_list, range_val=50):
    query = """
    with data(id) as (
        values %s
    )
    select s.id, s.mrdef, s.fda_labels, s.cas_reg_no, s.name, s.status
    FROM data d
    JOIN structures s 
    on s.id = d.id
    """
    cursor = getCursor()
    for i in range(0, len(dcid_list), range_val):
        chunk = dcid_list[i:i + range_val]
        chunk_tuple_list = tuple([(x,) for x in chunk])
        execute_values(cursor, query, chunk_tuple_list)
        for x in cursor:
            yield x


# In[28]:


list(bulkDCIDQueryIter([123,2446]))


# In[33]:


list(bulkSynQuerySubstringIter(["sirol","irol"]))


# In[34]:


def bulkDCIDCheckOralQueryIter(dcid_list, range_val=50):
    query = """
    with data(id) as (
        values %s
    )
    SELECT d.id
    FROM data d
    where EXISTS (SELECT 1 FROM struct2obprod so 
    JOIN ob_product op ON so.prod_id=op.id 
    WHERE op.route='ORAL'
    AND so.struct_id=d.id)
    """
    cursor = getCursor()
    for i in range(0, len(dcid_list), range_val):
        chunk = dcid_list[i:i + range_val]
        chunk_tuple_list = tuple([(x,) for x in chunk])
        execute_values(cursor, query, chunk_tuple_list)
        for x in cursor:
            yield x


# In[39]:


def testDCIDOral():
    for x in bulkDCIDCheckOralQueryIter(list(range(5392,5392+50))):
        print(x[0])


# In[40]:


def bulkDCIDGetSynonyms(dcid_list, range_val=50):
    query = """
    with data(id) as (
        values %s
    )
    select d.id, array_agg(s.lname)
    FROM data d
    join synonyms s on d.id=s.id 
    group by d.id
    """
    cursor = getCursor()
    for i in range(0, len(dcid_list), range_val):
        chunk = dcid_list[i:i + range_val]
        chunk_tuple_list = tuple([(x,) for x in chunk])
        execute_values(cursor, query, chunk_tuple_list)
        for x in cursor:
            yield x


# In[42]:


def testDCIDSynonyms():
    for x in bulkDCIDGetSynonyms(list(range(5392,5392+5))):
        print(x)


# In[45]:


def bulkDCIDGetApproval(dcid_list,range_val=50):
    query = """
    with data(id) as (
        values %s
    )
    select d.id, array_agg(a."type")
    FROM data d
    join approval a on d.id=a.struct_id 
    and (a."type"='EMA' or a."type"='FDA' or a."type"='PMDA')
    group by d.id
    """
    cursor = getCursor()
    for i in range(0, len(dcid_list), range_val):
        chunk = dcid_list[i:i + range_val]
        chunk_tuple_list = tuple([(x,) for x in chunk])
        execute_values(cursor, query, chunk_tuple_list)
        for x in cursor:
            yield x


# In[48]:


def testDCIDApproval():
    for x in bulkDCIDGetApproval(list(range(5392,5392+20))):
        print(x)


# In[ ]:




