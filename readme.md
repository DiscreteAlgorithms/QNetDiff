
This tool is QNetDiff by Shota nose et al. ([paper link])
You can download QNetDiff in (https://github.com/DiscreteAlgorithms/QNetDiff). 

QNetDiff is Python module, and it execute "python SparCC.py" command internal,  so it requires environment that "python" command is available. 

# Dependencys
 pip library name (test version)
 - matplotlib (3.5.1)
 - numpy (1.21.2)
 - pandas (1.3.5)
 - python-louvain (0.16)
 - scipy (1.7.3)

# Input

## 1. count_table_file

First row contains whole sample names, and first column contains whole bacteria names. Each sample or bacteira name must not contain space or tabs, and each name must be separetad by a space or tabs.
Other elements are a fraction representing the count of that row's bacteria read in that column's sample.

```
    sample_A    sample_B    ...
bacteria_a  0.2 ...  
bacteria_b  0.15
.               .
.                   .
.                       .
```

## 2. sample_to_group_file
Leftside element represents sample name and rightside element represents group name which the leftside sample is categorized. The num of samples (equals to num of rows) must be the same as count_table_file.


```
sample_A    group_A    
sample_B    group_B
sample_C    group_A
.
.
.
```  

## 3. category_file
Leftside represents bacteria name and rightside represents one level higher category (family for genus, genus for species, etc)

This information is used for unification of similar bacteria to reduce "false correlation" (ref [paper link]). 

```
Abiotrophia	Aerococcaceae
Acanthopleuribacter	Acanthopleuribacteraceae
Acaricomes	Micrococcaceae
.
.
.
```

# Usage
Execute the following commands in a shell that above library are available.

## Pattern 1 (basic)
```
python QNetDiff.py [count_table_file_path] [sample_to_group_file_path] [category_file_path] [focus_group] 
```
correlation file name automatically set to "[group_name].cor"

Pattern 1 example for example CRC data: 
```
python QNetDiff.py input/genus_count_table_H_S0.tsv input/sample_to_category_H_S0.tsv input/category.tsv Stage_0
```

## Pattern 2 (assigning correlation file name)
```
python QNetDiff.py [count_table_file_path] [sample_to_group_file_path] [category_file_path] [focus_group] [correlation_file_name]
```
correlation file name (SparCC output) set to [correlation_file_name]_[group].cor

## Pattern 3 (prepare correlation file and skip SparCC)
```
python QNetDiff.py [count_table_file_path] [sample_to_group_file_path] [category_file_path] [focus_group] [correlation_file_name] SKIP_SPARCC
```
skip SparCC and exist file (whose name is [correlation_file_name]_[group].cor) as correlation file. Even if you use exist correlation file, you need to prepare 3 input files. 

Example of Pattern 3 for CRC data: 
```
python QNetDiff.py input/genus_count_table_H_S0.tsv input/sample_to_category_H_S0.tsv input/category.tsv Stage_0 CRC_corr SKIP_SPARCC
```
