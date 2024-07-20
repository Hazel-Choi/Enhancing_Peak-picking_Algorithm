# Read the sample names from the text file
with open('unique_sample_names.txt','r') as file:
	sample_names=[line.strip() for line in file.readlines()]

# Create the download script
with open('download_all.sh', 'w') as script:
	script.write('#!/bin/bash\n\n')
	base_url = "https://github.com/Hazel-Choi/Enhancing_Peak-picking_Algorithm/releases/download/mzXML/"

	for sample in sample_names:
		for rep in ['A','B','C','D','E']:
			file_name=f"{sample}_{rep}.mzXML"
			download_url=f"{base_url}{file_name}"
			script.write(f'wget "{download_url}"\n')

print("Download script created: download_all.sh")