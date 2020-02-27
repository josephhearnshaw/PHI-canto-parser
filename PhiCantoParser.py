import json
import pandas as pd
from pandas.io.json import json_normalize
import re
from collections import OrderedDict


def metagenotypeParser(nodeData, curatedSession_Key, df):
    # Get the metagenotypes
    metagenotypes = [mg for mg in json_normalize(nodeData[curatedSession_Key]['metagenotypes'])]
    metagenotypes = [mg.split('.')[0] for mg in metagenotypes] # Remove the column names appended to the metagenotype identifiers
    metagenotype_df = pd.DataFrame(columns=['host_genotype', 'pathogen_genotype', 'type'])
    mg_list = []
    for metageno in list(OrderedDict.fromkeys(df['metagenotype'])):
        try:
            if metageno in metagenotypes:
                mg_list.append(metageno)
                metagenotype_df = metagenotype_df.append(pd.DataFrame(json_normalize(nodeData[curatedSession_Key]['metagenotypes'][metageno]), columns = metagenotype_df.columns), ignore_index=True)
        except:
            pass
    # Add the metagenotypes in, and walla! We have our pathogen and host genotypes. Now we need our pathogen genotype to be elucidated
    metagenotype_df['metagenotype'] = mg_list
    return metagenotype_df

with open('/home/joseph/phicanto.json', 'r') as f:
    allData = json.load(f)

### NEEDS TO BE PUT INTO A LOOP FOR ALL CURATION SESSIONS WHICH CAN BE INTERATED BASED ON THE LENGTH AND RANGE GIVEN BY NUMBER OF UNIQUE CURATION SESSIONS THEN DATAFRAMES MUST BE UNIQUELY NAMED
# Have all the curation sessions now
curationSessions = json_normalize(allData, ['curation_sessions'])
curatedSession_Keys = curationSessions[0]


# Now we need to get the value and normalize it
nodeData = {}
for key, node in allData['curation_sessions'].items():
    nodeData[key] = node



### NEEDS TO BE LOOPED FOR EACH CURATION SESSION AND CONCAT AT THE END
#Now have a loop that'll go over all the curation sessions, and create a datatable and write it out accordingly, we'll later concat this data into one big table
allele_data = pd.DataFrame(json_normalize(nodeData[curatedSession_Keys[0]]['alleles'], meta=['allele_type', 'gene', 'name', 'primary_identifier', 'synonyms', 'annotations'], ))
#
#metagenotypes = pd.DataFrame(json_normalize(nodeData[curatedSession_Keys[0]]['metagenotypes'], meta=['host_genotype', 'pathogen_genotype'] ))
metagenotypes = pd.DataFrame(json_normalize(nodeData[curatedSession_Keys[0]]['metagenotypes'] ))

### ANNOTATIONS
#Next we can get the annotations
annotations_df = pd.DataFrame(json_normalize(nodeData[curatedSession_Keys[0]], ['annotations']))
annotations_df = annotations_df.fillna('Not recorded')
annotations_df['Curation ID'] = curatedSession_Keys[0]
annotations_df['extension'] = annotations_df['extension'].apply(str).replace('[]', 'Not recorded')
# This will be the annotations dataframe with the gene ID recorded
annotations_df_nonExtensions = annotations_df[annotations_df.extension.str.contains('Not recorded') & (annotations_df.gene.str.contains('Not recorded') == False)]
# This is likely GOs related to the gene in question
annotations_df_nonExtensions = annotations_df_nonExtensions.drop(['conditions', 'metagenotype', 'genotype', 'extension'], 1)
# Separate the organism and gene names
annotations_df_nonExtensions['Organism'] = annotations_df_nonExtensions['gene'].str.rsplit().str[:-1]
annotations_df_nonExtensions['Organism'] = annotations_df_nonExtensions['Organism'].astype(str).str.replace(',', '').str.replace('[','').str.replace(']', '').str.replace('\'', '')
annotations_df_nonExtensions['gene'] = annotations_df_nonExtensions['gene'].str.split().str[-1]
cols_annotations = ['Curation ID', 'Organism', 'gene', 'type', 'term', 'term_suggestion.definition', 'term_suggestion.name', 'publication', 'evidence_code', 'submitter_comment', 'creation_date', 'curator.community_curated', 'curator.name', 'curator.email', 'checked' ]
annotations_df_nonExtensions = annotations_df_nonExtensions[cols_annotations]


#temporary - needs changing as it'll include multiple curation sessions
#annotations_df_nonExtensions.to_csv(f"/home/joseph/annotations_{curatedSession_Keys[0]}.txt", sep='\t', encoding='utf-8')



annotations_df_extensions_main = annotations_df[~annotations_df.extension.str.contains('Not recorded')]
annotations_df_extensions_main = annotations_df_extensions_main.reset_index()
# The extensions column doesn't use double qoutations and is therefore NOT a JSON.
annotations_df_extensions_main.extension = annotations_df_extensions_main.extension.str.replace('\'', "\"")

##### EXTENSIONS DATA FRAME #########
extensions_df = pd.DataFrame()
#del extensions_df
length, gene_list, term_list, mg_list = [], [], [], []
for index in range(len(annotations_df_extensions_main)):
    extensions_df = extensions_df.append(json_normalize(json.loads(annotations_df_extensions_main.extension[index])), sort=False)
    length.append(len(extensions_df))
    if index is not 0: # Appears to be the case that the length is 2 each time, but can't ascertain this will always be the case
        current_length = length[index] - length[index-1]
    else:
        current_length = length[index]
    [gene_list.append(annotations_df_extensions_main['gene'][index]) for _ in range(current_length)] # Append to the lists, to add later to the DF
    [term_list.append(annotations_df_extensions_main['term'][index]) for _ in range(current_length)]
    [mg_list.append(annotations_df_extensions_main['metagenotype'][index]) for _ in range(current_length)]

# Now we can add the respective terms and genes corresponding to their relative extensions
extensions_df['gene'], extensions_df['term'], extensions_df['metagenotype'], extensions_df['Curation ID'] = gene_list, term_list, mg_list, curatedSession_Keys[0]
extensions_df.fillna('Not recorded')

# Now we can write out the extensions, non extensions, pathogen host

#Obtain just hte pathogen host interaction
pathogen_host_interaction_df = annotations_df[annotations_df.type.str.contains("pathogen_host_interaction_phenotype") | annotations_df.type.str.contains("pathogen_phenotype") ]
pathogen_host_interaction_df.drop(['gene', 'genotype'], 1)

# Obtaining metagenotype data
metagenotype_df = metagenotypeParser(nodeData = nodeData, curatedSession_Key = curatedSession_Keys[0], df=pathogen_host_interaction_df)

pathogen_genotypes = metagenotype_df['pathogen_genotype']

genotype_df = pd.DataFrame(columns=['comment', 'loci', 'organism_strain', 'organism_taxonid'])
for i, pGenotype in enumerate(pathogen_genotypes):
    try:
        genotype_df = genotype_df.append(pd.DataFrame(json_normalize(nodeData[curatedSession_Keys[0]]['genotypes'][pGenotype]), columns = genotype_df.columns), ignore_index=True)
    except:
        pass

genotype_df['pathogen_genotype'] = pathogen_genotypes

# Obtaining the genotype mapping data for pathogens and hosts - this will be written out as a separate file in one big concat file with mapping data

###### METAGENOTYPE PATHOGEN HOST GENOTYPE MAPPING DATA ######
# Store list of list for pathogen IDs
pathogen_genoID = []
for i in genotype_df.loci:
    pathogen_genoID.append(i[0][0])

genotype_df['pathogen_geneID'] = pd.DataFrame(json_normalize(pathogen_genoID)['id'])
genotype_df = genotype_df.drop(['loci', 'pathogen_genotype'], 1)
genotype_df = genotype_df.fillna('Not recorded')
genotype_df.pathogen_geneID = genotype_df.pathogen_geneID.str.rsplit(':').str[0] # Gene name only
genotype_df = genotype_df[['pathogen_geneID', 'organism_strain', 'organism_taxonid', 'comment']] # recorder
genotype_df['metagenotype'] = metagenotype_df.metagenotype
filtered_metagenotype_df = genotype_df.merge(metagenotype_df, on='metagenotype') # Need metagenotype
filtered_metagenotype_df = filtered_metagenotype_df[['metagenotype', 'pathogen_genotype', 'pathogen_geneID', 'organism_strain', 'organism_taxonid', 'comment', 'host_genotype']] # Can individually write this out as mapping data

filtered_metagenotype_df

#### ALLELE DATA ####

alleles = pd.DataFrame(json_normalize(nodeData[curatedSession_Keys[0]], ['alleles']))

alleles = alleles.astype(str)
alleles = alleles[0].str.split(":", 1)

# Alleles table
allels_df = pd.DataFrame(nodeData[curatedSession_Keys[0]]['alleles'])
alleles_df_tf = allels_df.transpose()
alleles_df_tf['Curation Session'] = [curatedSession_Keys[0]][0] # CONSIDERATION - may need to map to allele name based on genotype


### MAY MAKE MORE SENSE TO PUT INTO ONDEX ---
## Or Parse thru and then get the annotations, add them. Then add the metagenotypes with Patho-Host interactions and types
#Get each allele
alleles = [x[0] for x in alleles_str]
