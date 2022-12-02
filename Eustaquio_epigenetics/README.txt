prepare your project folder:
1) create Procfile
    - "web: gunicorn app1:server"
2) create requirements.txt
    - this should include all libraries that your app needs to run on Heroku server

Before running app1 on orthofinder dashboard:
if copying old project directory, de sure to delete:
    1)  "./data/Orthofinder_data/summary_data/AssemblyAccession_to_SpeciesName.json"
    2)  "./data/all_orthogroup_info.json"
be sure that the file (1)AssemblyAccession_to_SpeciesName.json has beem made for the current dashboardspecies: 
Before uploading to Heroku servers:
1) make a virtual environment in your project folder and activate it
    $ python3 -m venv env
    $ source env/bin/activate
2) initialize a git repository
    $ git init
3) install all dependencies:
    $ pip3 install -r requirements.txt



```
app_folder
   |-- README.txt
   |-- assets
        |--style.css
        |--img.png
   |-- data
        |-- genome_annotations 
        |-- N0_HOG_counts.tsv
        |-- ncbi_dataset
            |-- data
                |--assembly_files...
        |-- New_proteomes 
        |-- Orthogroups
            |-- Orthogroups.GeneCount.tsv
            |-- Orthogroups.tsv 
        |-- Phylogenetic_Hierarchical_Orthogroups
        |-- Proteomes
            |-- proteome_files.faa...
            |-- Orthofinder
                |-- Results_MnthDay
                    |-- ...
        |-- Resolved_Gene_Trees
        |-- Species_Tree
        |-- summary_data
            |-- AssemblyAccession_to_SpeciesName.json
        ###the orthofinder dir should be added to the .gitignore file

        
        |-- summary_data
            |-- AssemblyAccession_to_SpeciesName.json

   |-- env
   |-- app1.py
   |-- Procfile
   |-- requirements.txt
   |-- .git
   |-- .gitignore

``` 

Steps for adding proteomes to orthofinder analysis:
1) get ncbi dataset proteomes (e.g. with proteomes_for_orthofinder.py)
2) create new directory in app_folder/data called 'New_proteomes'
    -copy proteomes from newly downloaded ncbi files into this folder and 
    rename with assembly accession #e.g. proteomes_for_orthofinder.grab_proteomes() will
    perform these steps
3) move ncbi assembly folder(s) with other older folders
4) run orthofinder
    $ orthofinder -f ./New_proteomes -b path/to/old/orthofinder/results_date
    
####Dash app preprocessing:
5) move all old folders and files produced from old runs to a folder named orthofinder_results_date
e.g. orthofinder_Results_old
    
    
#########need to update
5) add old results folder date to end of: 
    ./data/Orthofinder_data/summary_data/AssemblyAccession_to_SpeciesName.json, like below
    all_orthogroup_info.json      -->         all_orthogroup_info_Mar23_1.json
    AssemblyAccession_to_SpeciesName.json --> AssemblyAccession_to_SpeciesName_Mar23_1.json
6) rebuild data/Orthofinder_data/summary_data/AssemblyAccession_to_SpeciesName.json 

    