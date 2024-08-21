#!/usr/bin/env python
# coding: utf-8

import json5
import pandas as pd
from argparse import ArgumentParser

args = ArgumentParser('./build_bwa_scripts.py', description="""This program has been designed to
    generate bash scripts to align genome sequencing data using BWA, Samtools, and PicardTools. See
    documentation for links to the respective tools. """)

args.add_argument(
    '-c',
    '--config_file',
    help="""
    This is the config json file containing information about parameters and links to necessary programs
    """,
    default="bwa_parameters.json"
)

args.add_argument(
    '-s',
    '--samples',
    help="""
    This is the file containing information about sample names and fastq file information""",
    default="samples.txt"
)

args = args.parse_args()

with open(args.config_file, 'r') as cf:
    config_parameters = json5.load(cf)

samples_df = pd.read_csv(args.samples)

def build(df, config):
    output_info = {}
    for index, row in df.iterrows():
        output_commands = []
        num_cores = str(config_parameters["number_cores"])
        fastq_variable = row[1]
        if config['is_paired']:
            fastq_variable=row[1] +' '+row[2]
        output_commands.append('## Running on '+num_cores+' threads')
        output_commands.append('## The first step is to align to the genome')
        # The first step is to align the file. How this command is formatted will depend on several parameters
        if config['read_length'] > 69:
            # read length suggests that using bwa mem will be more favorable
            output_commands.append('## Read length suggests that bwa mem is the more favorable option')
            output_commands.append(config['path_to_bwa']+' mem -t '+num_cores+' '+config['path_to_prebuilt_reference']+' '+fastq_variable+' > '+row[0]+'.sam')
        output_commands.append('## Checking alignment stats')
        output_commands.append(config['path_to_samtools']+' stats '+row[0]+'.sam'+' | grep ^SN | cut -f 2- > '+row[0]+'.01sam.align_stats.txt')
        output_commands.append('## Save unmapped reads and convert them to fastq format for troubleshooting')
        output_commands.append(config['path_to_samtools']+' view -b -f 4 -@ '+num_cores+' -t '+config['fai_file']+' '+row[0]+'.sam > '+row[0]+'.unmapped.bam')
        output_commands.append(config['path_to_samtools']+' fastq '+row[0]+'.unmapped.bam -1 '+row[0]+'_unmapped1.fastq -2 '+row[0]+'_unmapped2.fastq -s '+row[0]+'_single_read_mapped.fastq')
        output_commands.append('## Convert to BAM format and remove multimapped reads')
        output_commands.append(config['path_to_samtools']+' view -bq 10 -@ '+num_cores+' '+row[0]+'.sam'+' > '+row[0]+'.uniquelyAligned.bam')
        output_commands.append('## Remove large SAM file and calculate stats for uniquely aligned')
        output_commands.append('rm '+row[0]+'.sam')
        output_commands.append(config['path_to_samtools']+' stats '+row[0]+'.uniquelyAligned.bam | grep ^SN | cut -f 2- > '+row[0]+'.02uniquelyAligned.align_stats.txt')
        output_commands.append('## Now to label the sample inside the bam file, required for Picard Tools')
        output_commands.append(config['path_to_samtools']+' addreplacerg -@ '+num_cores+' -r "@RG\\tID:RG1\\tSM:'+row[0]+'\\tPL:Illumina\\tLB:Library.fa" -o '+row[0]+'.uniquelyAlignedLabeled.bam '+row[0]+'.uniquelyAligned.bam')
        output_commands.append('## Now sort it and index for Picard Tools')
        output_commands.append(config['path_to_samtools']+' sort -@ '+num_cores+' '+row[0]+'.uniquelyAlignedLabeled.bam -o '+row[0]+'.uniquelyAlignedLabeledSorted.bam')
        output_commands.append(config['path_to_samtools']+' index -@ '+num_cores+' '+row[0]+'.uniquelyAlignedLabeledSorted.bam')
        output_commands.append('## Now use PicardTools to mark and remove duplicates')
        output_commands.append('java -jar '+config['path_to_picard_jar']+' MarkDuplicates -REMOVE_DUPLICATES true -VALIDATION_STRINGENCY SILENT -AS true -I '+row[0]+'.uniquelyAlignedLabeledSorted.bam -O '+row[0]+'.rmdup.bam -M '+row[0]+'.picardMetrics.txt')
        output_commands.append('## Run final stats and index the bam file')
        output_commands.append(config['path_to_samtools']+' index -@ '+num_cores+' '+row[0]+'.rmdup.bam')
        output_commands.append(config['path_to_samtools']+' stats '+row[0]+'.rmdup.bam'+' | grep ^SN | cut -f 2- > '+row[0]+'.03rmdup.align_stats.txt')
        # Now save everything to the output directory
        output_info[row[0]] = output_commands
    return output_info

output_parameters = build(samples_df, config_parameters)
for sample in output_parameters:
    with open(sample+'_bwa.sh', 'w') as script_file:
        for item in output_parameters[sample]:
            script_file.write(f"{item}\n")
