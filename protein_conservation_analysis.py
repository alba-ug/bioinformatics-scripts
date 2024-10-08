#!/usr/bin/python3.8
import os, sys, shutil, subprocess, re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#START of Funtions:

#Function to validate working directory
def dir_input() :
  while True :
    print("\nWARNING: If a directory with the same name is present in the given path, the program will delete it.\n")
    user_dir = input("Enter the path and the name of the directory to store the program files:\n(The new directory name must contain only letters, numbers or _, e.g., /home/ICA2/Protein_analys)\n\t")
    initialwrkdir = os.getcwd()
    parent_cdir = re.match('/.+/',initialwrkdir)
    local_user_dir = "/localdisk" + user_dir
#Prevent user from deleting current working directory
    if initialwrkdir == user_dir or initialwrkdir == local_user_dir :
      print("\nThe input is not valid. The current program would be deleted.")
#Prevent user from deleting parent directories from current directory
    elif re.search(user_dir, parent_cdir[0][:-1]) :
      print("\nThe input is not valid. The current program would be deleted.")
    elif os.path.exists(user_dir) and initialwrkdir != user_dir and initialwrkdir != local_user_dir :
      shutil.rmtree(user_dir)
      os.mkdir(user_dir)
      os.chdir(user_dir)
      break
    else :
#The new directory must contain only letters, numbers and _ 
      if re.fullmatch('(\/.+)(\/[a-zA-Z0-9_]+)+',user_dir) :
         try :
           parent_dir = re.fullmatch('(\/.+)(\/[a-zA-Z0-9_]+)+',user_dir)[1]
           os.chdir(parent_dir)
           new_dir = re.fullmatch('(\/.+)(\/[a-zA-Z0-9_]+)+',user_dir)[2]
           new_dir = new_dir[1:]
           os.mkdir(new_dir)
           os.chdir(new_dir)
           user_dir = os.getcwd()
           break
         except :
           print("\nThe input is not valid.")
      else :
         print("\nThe input is not valid.")
  return user_dir


#Functions to evaluate user input for taxon name or taxon ID.
#It accepts a single word(can contain - or . if more words are present) or words separated by one space to specify the subset of thetaxonomic tree.
def usertaxon() :
	utaxon = input ("\nWhich subset of the taxonomic tree do you want to analyze?\n(e.g., birds, or g-proteobacteria)\n\t")
	while re.fullmatch('[a-zA-Z]+(\-[a-zA-Z])*', utaxon) == None and re.fullmatch('[a-zA-Z]+\.?(\s\S+\.?\S+)+', utaxon) == None:
		utaxon = input ("Please enter a valid taxon name:\n\t")
	utaxon = utaxon.lower().replace(" ", "_")
	return utaxon

#It accepts non negative integer
def usertaxid() :
        userid = input ("\nEnter the ID of the taxonomic subset you would like to analyze:\n\t")
        while re.fullmatch('[0-9]+', userid) == None:
                userid = input ("Please enter a valid taxon id:\n\t")
        return userid

#Functions to evaluate query results, first one looks for taxon name and second one taxon id, they include error trap for no results and NCBI query error
def esearch_namerror() :
	tnames = ""
#Ask user to provide taxonomic group until successful query
	while tnames == "" or "WARNING" in tnames or "ERROR" in tnames or "FAILURE" in tnames :
		utaxon = usertaxon()
		print("\nSearching for taxon name matching input...")
		query_taxname = "esearch -db taxonomy -spell -query " + utaxon + " |efetch "
		tnames = subprocess.getoutput(query_taxname)
#In case bad functionaly from NCBI website
		while "Bad Request" in tnames :
			tnames = subprocess.getoutput(query_taxname)
		if tnames == "" or "WARNING" in tnames or "ERROR" in tnames :
			print("\nSorry, there's no associated taxon with your input.")
	query_taxid = "esearch -db taxonomy -query " + utaxon + " |efetch -format taxid "
	tids = subprocess.getoutput(query_taxid)
	return tnames, tids


def esearch_iderror() :
	tnames = ""
#Ask user to provide NCBI taxon ID until successful query
	while tnames == "" or "WARNING" in tnames or "ERROR" in tnames or "FAILURE" in tnames :
                userid = usertaxid()
                print("\nSearching for taxon ID matching input...")
                query_taxid = "esearch -db taxonomy -query " + userid + "[uid] |efetch "
                tnames = subprocess.getoutput(query_taxid)
#In case bad functionaly from NCBI website
                while "Bad Request" in tnames :
                        tnames = subprocess.getoutput(query_taxid)
                if tnames == "" or "WARNING" in tnames or "ERROR" in tnames :
                        print("\nSorry, there's no associated taxon with your input.")
	tids = userid
	return tnames, tids


#Function to validate integer input [0, inf)
def int_input() :
  while True :
    numinput = input("\nEnter a number from the given options:\n\t")
    if numinput.isdigit() :
        break
    else :
        print("The answer is not valid.")
  numinput = int(numinput)
  return numinput

#Function to get a number from the user in a specific range
def select_num(definedlist) :
  unumber = int_input()
  while unumber not in definedlist :
    unumber = int_input()
  return unumber


#The function asks the user to input taxon name. The function searches on NCBI taxonomy for the user input and the user selects from the given options.
#This process will repeat itself until the user agrees with the search results.
def get_taxon() :
  confirmation = ""
  taxname = ""
  taxid = ""
  while confirmation != "y":
    print("\nWould you like to provide the taxon name or the NCBI taxon id?\n\t1. Taxon name\n\t2. NCBI taxon ID")
    initial_options = [1,2]
    preference = select_num(initial_options)
    if preference == 1 :
#Call usertaxon function to get user input inside esearch_namerror to get a usable taxon name
      tnames,tids = esearch_namerror()
    elif preference == 2 :
#Call usertaxid function to get user input inside esearch_iderror to get a usable taxon ID
      tnames,tids = esearch_iderror()
#Evaluates if there's more than one option and lets the user choose
    if tids.count("\n") > 0 :
      max = int((tnames.count("\n")+1)/2)
      print("\nWhich option would you like to pursue?\n\t")
      print(tnames)
      optionslist = list(range(1, max+1))
      userselection = select_num(optionslist)
      taxid = tids.split("\n")[userselection-1]
      i1 = (2*userselection)-2
      taxname = tnames.split("\n")[i1].split(".")[1]
      taxname = taxname.replace(" ","_")
      comname = tnames.split("\n")[i1+1].split(",")[0]
      comname = comname.replace(" ","_")
      decision = input ("\nDo you want to analyze " + taxname[1:] + comname[4:] + "? y/n\n\t")
      confirmation = decision.lower().replace(" ","")
#When there's only one result in the search, asks the user to confirm in order to proceed
    else :
      taxid = tids
      taxname = tnames.split()[1]
      taxname = taxname.replace(" ","_")
      comname = tnames.split()[2]
      comname = comname.replace(" ","_")
      decision = input ("\nDo you want to analyze " + taxname + comname[:-1] + "? y/n\n\t")
      confirmation = decision.lower().replace(" ","")
  return taxname, taxid


#Function to retrieve protein family name, it accepts spaces, -, and alphanumeric characters
def get_protname() :
        protname = input ("\nEnter the name of the protein family you would like to analyze:\n(write the protein family name in the singular form, e.g., catalase)\n\t")
        while re.fullmatch('[a-zA-Z]+[ \w\-]*', protname) == None :
                protname = input ("Please enter a valid protein family name:\n\t")
        return protname


#Function to determine the number of protein sequences with trap error for failure in NCBI query
def num_seqs(protname,taxid) :
    print("\nSearching for sequences to analyze...\nWait for confirmation to proceed with analysis.")
    query = 'esearch -db protein -spell -query \"' + protname + ' [PROT] AND txid' + taxid + ' [ORGN] NOT (predicted OR hypothetical)\"'
    textproc = query + ' | grep "Count"'
    search_count = subprocess.getoutput(textproc)
    while "FAILURE" in search_count :
       search_count = subprocess.getoutput(textproc)
    total_seqs = int(search_count.replace(" ","").replace("<Count>","").replace("</Count>",""))
    return total_seqs


#Function to ask user how to proceed with analysis [Continue, Modify query, Exit]
def get_confirmation() :
      print("\nHow do you want to proceed?\n\t1. Continue with analysis.\n\t2. Modify query.\n\t3. Exit program.")
      initial_options = [1,2,3]
      preference = select_num(initial_options)
      if preference == 1 :
        confirmation = "y"
      elif preference == 2 :
        confirmation = "n"
      else :
        exit()
      return confirmation


#Function to get taxonid, protein family name and user agreement to proceed with number of sequences
def user_input() :
  confirmation = ""
  while confirmation != "y":
#Ask user for taxon name or taxon ID and protein family
    taxname, taxid = get_taxon()
    protname = get_protname()
    basename = protname.replace(" ","_") + "_" + taxname.replace(" ","_") + "_taxID" + taxid
#Find the number of protein sequences given the search
    total_seqs = num_seqs(protname,taxid)
#Evaluates the number of sequences to promt message for user and asks confirmation to proceed
    if total_seqs == 0 or total_seqs == 1 :
      print("\nNo results were found for the query.")
      print("\nHow do you want to proceed?\n\t1. Modify query.\n\t2. Exit program.")
      initial_options = [1,2]
      preference = select_num(initial_options)
      if preference == 2 :
        exit()
    elif total_seqs > 1 and total_seqs < 5 :
      print("\nWARNING: There are less than 5 sequences associated with the query.")
      decision = input ("\nDo you want to proceed with the analysis? y/n\n\t")
      confirmation = decision.lower().replace(" ","")
    elif total_seqs > 1000 :
      print("\nWARNING: There are more than 1000 sequences associated with the query.")
      print("The program can analyze up to 1000 sequences.")
      print("Refine your search or continue analysis with the first 1000 sequences in the query.")
      confirmation = get_confirmation()
    else :
      print("\nThere are " + str(total_seqs) + " sequences associated with the query.")
      confirmation = get_confirmation()
  return taxname, taxid, protname, basename, total_seqs


#Function to get protein sequences given taxid and protein family, max 1000 sequences
def get_protseq(taxid,protname,basename,total_seqs) :
    filename = basename + ".fasta"
    query = 'esearch -db protein -query \"' + protname + ' [PROT] AND txid' + taxid + ' [ORGN] NOT (predicted OR hypothetical)\"'
    if total_seqs > 1000 :
      subptext = query + " | efetch -format fasta -stop 1000 "
    else :
      subptext = query + " | efetch -format fasta "
    fastasequences = subprocess.getoutput(subptext)
    while "FAILURE" in fastasequences :
      fastasequences = subprocess.getoutput(subptext)
    with open(filename,"w") as file:
     file.write(fastasequences)
    return fastasequences

#Function to get protein sequences given taxid and protein family, max 1000 sequences
def get_accessions(taxid,protname,total_seqs) :
    query = 'esearch -db protein -query \"' + protname + ' [PROT] AND txid' + taxid + ' [ORGN] NOT (predicted OR hypothetical)\"'
    if total_seqs > 1000 :
      subptext = query + " | efetch -format docsum -stop 1000 | xtract -pattern DocumentSummary -element OSLT"
    else :
      subptext = query + " | efetch -format docsum | xtract -pattern DocumentSummary -element OSLT "
    accessions = subprocess.getoutput(subptext)
    while "FAILURE" in accessions :
      accessions = subprocess.getoutput(subptext)
    return accessions

#Function to get protein sequences given taxid and protein family, max 1000 sequences
def get_species(taxid,protname,total_seqs) :
    query = 'esearch -db protein -query \"' + protname + ' [PROT] AND txid' + taxid + ' [ORGN] NOT (predicted OR hypothetical)\"'
    if total_seqs > 1000 :
      subptext = query + " | efetch -format docsum -stop 1000 | xtract -pattern DocumentSummary -element Organism "
    else :
      subptext = query + " | efetch -format docsum | xtract -pattern DocumentSummary -element Organism "
    species = subprocess.getoutput(subptext)
    while "FAILURE" in species :
      species = subprocess.getoutput(subptext)
    return species

# END of Funtions


#INTRODUCTION TO PROGRAM

print("\n------------------------------------------------------------------------------------------------------------\n")
print("Welcome,\n\nThe purpose of this program is to perform a protein sequence conservation analysis.\n")
print("You will be asked to provide the name of a subset of the taxonomic tree and the name of a protein family you are interested in.\n")
print("The program will also ask to provide a path to a directory to store the results.\n")
print("The program will carry out the following analysis:\n\n")
print("\t- Protein sequence conservation analysis (EMBOSS program: plotcon)\n")
print("\t- Identification of motifs from the PROSITE database (EMBOSS program: patmatmotifs)\n")
print("\t- Prediction of protein secondary structure (GOR method) (EMBOSS program: garnier)\n")
print("\t- Phylogenetic tree construction (PhyML)\n")
print("\n------------------------------------------------------------------------------------------------------------\n")


# ICA Point 1 - Get user input

#Create directory to store program files, user decides path and name or done automatically
print("\nA directory will be created to store all the files required for the program and to store the results.")
print("\nHow would you like to proceed?\n\t1. The program will create the Protein_Analysis directory in the current directory\n\t2. Enter the path and name of the new directory")
print("\nWARNING: If a directory with the same name is present in the given path, the program will delete it.")
initial_options = [1,2]
preference = select_num(initial_options)
if preference == 1 :
   wrkdir = os.getcwd()
   os.chdir(wrkdir)
   try:
#If there is a directory in the chosen path with the name Protein_analysis, it will remove it.
      shutil.rmtree("Protein_analysis")
      os.mkdir("Protein_analysis")
      os.chdir("Protein_analysis")
   except :
      os.mkdir("Protein_analysis")
      os.chdir("Protein_analysis")
   wrkdir = os.getcwd()
else :
   wrkdir = dir_input()


#Ask user if agrees with number of species and sequences to be analyzed given a taxon and protein family
final_confirmation = ""
while final_confirmation != "y":
#Run functios to get all user input
    taxname,taxid,protname,basename,total_seqs= user_input()
#Get protein sequences given user input, accession numbers and species
    print("\nWait for confirmation to proceed with analysis.")
    species = get_species(taxid,protname,total_seqs)
    list_species = species.split("\n")
    total_species = len(set(list_species))
#Evaluates the number of sequences to promt message for user and asks confirmation to proceed
    if total_species < 5 :
      print("\nWARNING: There are less than 5 species in the dataset.")
      final_confirmation = get_confirmation()
    elif total_species > 200 :
      print("\nWARNING: There are " + str(total_species) + " species in the dataset.")
      final_confirmation = get_confirmation()
    else :
      print("\nThere are " + str(total_seqs) + " sequences from " + str(total_species) + " species in the dataset.")
      final_confirmation = get_confirmation()

#Ask user if they want to include simple post-translational modifications in the motifs scanning
print("\nWhile determining whether any known motifs are associated with the sequences, do you want to include simple post-translational modifications in the results?\n\t1. Yes\n\t2. No")
initial_options = [1,2]
prune = select_num(initial_options)

print("\nSTARTING ANALYSIS...\n")
print("Depending on the size of the dataset, it can take up to several minutes to complete all analysis.\nPreliminary results will be shown on the screen.")
print("\n------------------------------------------------------------------------------------------------------------\n")

#Once user input is complete, retrieve sequences and accession numbers
sequences = get_protseq(taxid,protname,basename,total_seqs)
filename = basename + ".fasta"
#sequences = open(filename).read()
accessions = get_accessions(taxid,protname,total_seqs)
list_accessions = accessions.split("\n")

# ICA Point 2

print("STARTING PROTEIN SEQUENCE CONSERVATION ANALYSIS\n")
#Perform sequence alignment and create conservation plot
subprocess.call("clustalo -i " + basename + ".fasta -o Protein_alignment.msf --outfmt msf --output-order=tree-order --residuenumber --threads=200 --force ", shell=True)
subprocess.call("plotcon -sequences Protein_alignment.msf -winsize 4 -graph pdf -goutfile " + basename + "_plotcon", shell=True)
subprocess.call("plotcon -sequences Protein_alignment.msf -winsize 4 -graph png -goutfile " + basename + "_plotcon", shell=True)
print("\nThe conservation plot will be displayed on screen once ALL analysis are FINISHED.")
print("\nPROTEIN SEQUENCE CONSERVATION ANALYSIS FINISHED")
print("\n------------------------------------------------------------------------------------------------------------\n")


# ICA Point 3

print("STARTING MOTIFS SCANNING\n")
os.mkdir("Patmatmotifs_result_files")
os.chdir("Patmatmotifs_result_files")
#Get individual headers
fasta_headers = re.findall('>.+\n',sequences)
#Get individual sequences
individual_seqs = re.split('>.+\n',sequences)
individual_seqs.pop(0)

list_seqnum = []
total_motifs = {}
for i in range(len(list_accessions)) :
 # print(i)
  sequencenum = "Seq_" + str(i+1)
  list_seqnum.append(sequencenum)
  filename1 = "seq" + str(i+1) + ".fasta"
  filename2 = list_accessions[i] + ".patmatmotifs"
#Create individual fasta files for each sequence
  with open(filename1,"w") as file:
     file.write(fasta_headers[i] + individual_seqs[i])
#Create .patmatmotifs files with or without pruning
  if prune == 1 :
     subprocess.call("patmatmotifs -sequence " + filename1 + " -outfile " + filename2 + " -prune N -auto", shell=True)
  elif prune == 2 :
     subprocess.call("patmatmotifs -sequence " + filename1 + " -outfile " + filename2 + " -auto", shell=True)
#Remove .fasta files
  os.remove(filename1)
  motifs_file = open(filename2).read()
  with open("Associated_motifs.txt","a") as file:
     file.write(motifs_file)
  motifs = re.findall('Motif = .+\n',motifs_file)
  motifs_list = []
  for motif in motifs :
    m = motif.replace("Motif = ","").replace("\n","")
    motifs_list.append(m)
#Create dictionary with sequence number and the motifs they contain
  total_motifs[sequencenum] = motifs_list
path = wrkdir + "/Associated_motifs.txt"
os.rename("Associated_motifs.txt", path)
os.chdir(wrkdir)

#Get the set of motifs considering all the sequences
unique_motifs = []
for motifs_list in total_motifs.values():
    for motif in motifs_list :
        unique_motifs.append(motif)
set_motifs = sorted(list(set(unique_motifs)))

print("The folowing motifs where found in the dataset:\n")
for motif in set_motifs :
   print("\n\t- " + motif)

#Create dictionary with sequence number and the motifs counts
counts_motifs = {}
for sequencenum, motifs_list in total_motifs.items() :
    counts_list = []
    for motif in set_motifs :
        count = motifs_list.count(motif)
        counts_list.append(count)
    counts_motifs[sequencenum] = counts_list

#Convert dictionary with motifs counts into dataframe
df = pd.DataFrame()
df1 = df.from_dict(counts_motifs)
df_motifs = pd.DataFrame()
df_motifs = df1.transpose()
#Rename dataframe columns to mofif names
for i in range(len(set_motifs)) :
   df_motifs = df_motifs.rename(columns={ i : set_motifs[i] })

#Graph with total counts per motif
motifs_x = list(df_motifs.columns)
plt.bar(motifs_x,df_motifs.sum(),color='midnightblue')
plt.title("Motifs scanning (PROSITE)")
plt.xlabel('Motif')
plt.ylabel('Total count per motif')
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.savefig('Motifs_counts.png')
plt.close()

#Add sequence identifiers to data frame
final_df_motifs = df_motifs
#Create dictionaries with sequence number as key and species and accessions as values
dic_seqnum_species = {}
for i in range(len(list_seqnum)) :
   dic_seqnum_species[list_seqnum[i]] = list_species[i]

dic_seqnum_accs = {}
for i in range(len(list_seqnum)) :
   dic_seqnum_accs[list_seqnum[i]] = list_accessions[i]

final_df_motifs['Accession'] = pd.Series(dic_seqnum_accs)
final_df_motifs['Species'] = pd.Series(dic_seqnum_species)
#final_df_motifs.to_csv("Motifs_counts.tsv",sep="\t",header=True)
final_df_motifs.to_csv("Motifs_counts.csv",header=True)

print("\nSummary of motifs found in dataset:\n")
print(final_df_motifs.describe())

if total_seqs > 100 :
  final_df_motifs = final_df_motifs.iloc[0:100,:]

#Graph with counts per motif per sequence
os.mkdir("Motifs_counts_graphs")
os.chdir("Motifs_counts_graphs")
for motif in motifs_x :
  plt.figure(figsize=(20,10))
  plt.bar(final_df_motifs['Accession'], final_df_motifs[motif], label = motif, width=0.2)
  name_motif = motif.replace(" ","_")
  plt.title('Counts of ' + motif + ' per sequence')
  plt.xlabel('Accession number')
  plt.ylabel('Counts')
  plt.xticks(rotation=90, ha='right')
  plt.tight_layout()
  plt.legend()
  plt.savefig(name_motif + '_counts_.png')
  plt.close()
os.chdir(wrkdir)

print("\nThe following files were created:\n\n\t- Patmatmotifs_result_files (directory with files with motif scanning results for each sequence)\n")
print("\t- Associated_motifs.txt (file with motif scanning with results from all sequences)\n")
print("\t- Motifs_counts.csv (file with data frame with counts per motif per sequence)\n")
print("\t- Motifs_counts.png (file with plot of total motifs counts)\n")
print("\t- Motifs_counts_graphs (directory with plots of motifs counts per sequence, up to 100 sequences)\n")
print("\nMOTIFS SCANNING FINISHED")
print("\n------------------------------------------------------------------------------------------------------------\n")


# ICA Point 4

print("STARTING SECONDARY STRUCTURE PREDICTION\n")
subprocess.call("garnier -sequence " + basename + ".fasta -outfile " + basename + "_secondary_structure.garnier", shell=True)
print("\nThe following file was created:\n\n\t- " + basename + "_secondary_structure.garnier (file with secondary structure prediction for each sequence)\n")
print("\nSECONDARY STRUCTURE PREDICTION FINISHED")
print("\n------------------------------------------------------------------------------------------------------------\n")

# ICA Point 4

print("\nSTARTING PHYLOGENETIC TREE ANALYSIS")
#Create names for sequences for phylogenetic tree
phy_names = []
i = 1
for sp in list_species:
 nodotsp = sp.replace(".","")
 try:
  first = nodotsp.split()[0]
  initial = list(first)[0]
  second = nodotsp.split()[1]
  if total_seqs < 100 :
    partialname = second[:5]
  elif total_seqs > 99 and total_seqs < 1000 :
    partialname = second[:4]
  elif total_seqs >= 1000 :
    partialname = second[:3]
  name = initial + "_" + partialname + str(i)
  phy_names.append(name)
  i += 1
 except :
  if len(list(nodotsp)) < 6 :
    partialname = nodotsp
  else :
    partialname = nodotsp[:5]
  name = partialname + str(i)
  phy_names.append(name)
  i += 1

#Create list of modified accessions (replace . by _)
modified_acc = []
for acc in list_accessions :
   new_acc = acc.replace(".","_")
   modified_acc.append(new_acc)
#Create list of accessions ordered by alignment
alignment = open("Protein_alignment.msf").read()
accs_align = re.findall('Name.+\n',alignment)
ordered_accs = []
for acc in accs_align :
   ordered_accs.append(acc.split()[1])
if total_seqs > 100 :
   ordered_accs_top = ordered_accs[:100]
else :
   ordered_accs_top = ordered_accs

#Select sequences according to clustal alignment
for i in range(len(modified_acc)) :
  if modified_acc[i] in ordered_accs_top:
#Create new fasta file for use in alignment for phyml
     name = ">" + phy_names[i] + "\n"
     with open("Topseqs_phy.fasta","a") as file:
        file.write(name + individual_seqs[i])

#Alingment with max 100 sequences
subprocess.call("clustalo -i Topseqs_phy.fasta -o Protein_alignment.phy --outfmt phy --output-order=tree-order --threads=200 --force ", shell=True)
#Tree creation with phyml
subprocess.call("phyml -i Protein_alignment.phy -d aa ", shell=True)
#subprocess.call("phyml -i Protein_alignment.phy -d aa > Tree_information.txt", shell=True)
subprocess.call("figtree -graphic PDF Protein_alignment.phy_phyml_tree.txt " + basename + "_Tree.pdf", shell=True)
subprocess.call("figtree -graphic PNG -width 800 -height 500 Protein_alignment.phy_phyml_tree.txt " + basename + "_Tree.png", shell=True)
print("\nThe following files were created:\n\n\t- Protein_alignment.phy (alignment of sequences in phy format)\n")
print("\t- Protein_alignment.phy_phyml_tree.txt (tree file)\n")
print("\t- Protein_alignment.phy_phyml_stats.txt (model parameters file)\n")
print("\t- " + basename + "_Tree.pdf (tree image file)\n")
print("\t- Tree_information.txt (information on tree creation process)\n")
print("\nPHYLOGENETIC TREE ANALYSIS FINISHED")
print("\n------------------------------------------------------------------------------------------------------------\n")

# ICA Point 2 - SHOW ON SCREEN
print("STARTING RESULTS DISPLAY\n")
print("It can take up to minutes for the plots to display, please wait.")
print("\n\t(Image 1/3) Displaying protein conservation plot ...")
subptext = "eog " + basename + "_plotcon.1.png"
subprocess.call(subptext, shell=True)
print("\n\t(Image 2/3) Displaying motifs found in dataset ...")
subprocess.call("eog Motifs_counts.png", shell=True)
print("\n\t(Image 3/3) Displaying phylogenetic tree ...")
subptext = "eog " + basename + "_Tree.png"
subprocess.call(subptext, shell=True)
print("\nRESULTS DISPLAY FINISHED")
print("\n------------------------------------------------------------------------------------------------------------\n")
print("EXITING PROGRAM\n")
os.remove("Topseqs_phy.fasta")
