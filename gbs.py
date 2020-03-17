#!/usr/local/env python3

# Dependencies
# conda install bbmap
# conda install pigz
# pip install pysam
# conda install blast
# conda install bowtie2
# conda install samtools
# conda install kmc
# conda install jellyfish
# conda install psutil
# conda install biopython
# conda install natsort
# conda install pyVCF


from gbs_methods import Methods
from argparse import ArgumentParser
from multiprocessing import cpu_count
from psutil import virtual_memory
import os
from glob import glob


class GBS(object):
    def __init__(self, args):
        # Command line arguments
        self.input = args.input
        self.out_folder = args.output
        self.cpu = args.threads
        self.mem = args.memory
        self.ref = args.reference
        if args.paired_end:
            self.read_type = 'pe'
        elif args.single_end:
            self.read_type = 'se'
        elif (args.paired_end and args.single_end)\
                or (not args.paired_end and not args.single_end):
            raise Exception('Please choose between "-se" and "-pe"\n')
        self.adapters = args.adapters

        # Data
        self.sample_dict = dict()

        # Run
        self.run()

    def run(self):
        print('Checking a few things...')
        # Check if number of CPU and memory requested are valid
        self.cpu = Methods.check_cpus(self.cpu)
        self.mem = Methods.check_mem(self.mem)

        # Check if folders are not empty
        result = Methods.check_folder_empty(self.input)
        if result == 0:
            raise Exception('Input folder does not contain files with accepted file extensions: {}'.format(
                Methods.accepted_extensions))

        # Get input files and place info in dictionary
        Methods.get_files(self.input, self.sample_dict)

        # Trim reads
        adapter_dict = dict()
        # if self.adapters:
        #     Methods.parse_adapters(self.adapters, adapter_dict)
        #     if any(info.direction == 'left' for ident, info in adapter_dict.items())\
        #             and not any(info.direction == 'right' for ident, info in adapter_dict.items()):
        #         Methods.parallel_trim_reads(Methods.trim_left, adapter_dict, self.sample_dict,
        #                                     self.out_folder + '/trimmed/', self.mem, self.cpu)
        #     elif not any(info.direction == 'left' for ident, info in adapter_dict.items())\
        #             and any(info.direction == 'right' for ident, info in adapter_dict.items()):
        #         Methods.parallel_trim_reads(Methods.trim_right, adapter_dict, self.sample_dict,
        #                                     self.out_folder + '/trimmed/', self.mem, self.cpu)
        #     elif any(info.direction == 'left' for ident, info in adapter_dict.items())\
        #             and any(info.direction == 'right' for ident, info in adapter_dict.items()):
        #         Methods.parallel_trim_reads(Methods.trim_both, adapter_dict, self.sample_dict,
        #                                     self.out_folder + '/trimmed/', self.mem, self.cpu)

        # Index reference genome if not already done
        ref_index = '.'.join(self.ref.split('.')[:-1])
        if not os.path.exists(ref_index + '.1.bt2'):
            Methods.index_bowtie2(self.ref, ref_index, self.cpu)

        # Map fastq on indexed reference genome
        Methods.parallel_map_bowtie2(self.out_folder + '/mapped/', ref_index, self.sample_dict, self.cpu)

        # Call variants
        if not os.path.exists(self.ref + '.fai'):
            Methods.index_samtools(self.ref)

        # List bam files into a text file, one per line
        bam_list = glob(self.out_folder + '/mapped/' + '*.bam')
        Methods.list_to_file(bam_list, self.out_folder + '/bam.list')

        print('Calling variants with bcftools...')
        Methods.call_variants(self.ref, self.out_folder + '/bam.list',
                              self.out_folder + '/variants.vcf', self.cpu)
        # Sort VCF file
        print('Fixing vcf file...')
        Methods.fix_vcf(self.out_folder + '/variants.vcf',  self.out_folder + '/variants_fixed.vcf')

        print('Sorting vcf file...')
        Methods.sort_vcf(self.out_folder + '/variants_fixed.vcf',
                         self.out_folder + '/variants_sorted.vcf')


if __name__ == "__main__":
    max_cpu = cpu_count()
    max_mem = int(virtual_memory().total * 0.85 / 1000000000)  # in GB

    parser = ArgumentParser(description='Reference-based genotyping-by-sequencing (GBS) pipeline')
    parser.add_argument('-r', '--reference', metavar='/reference_genome.fasta',
                        required=True,
                        help='Reference genome for read mapping')
    parser.add_argument('-i', '--input', metavar='/input_folder/',
                        required=True,
                        help='Folder that contains the fastq files')
    parser.add_argument('-o', '--output', metavar='/output_folder/',
                        required=True,
                        help='Folder to hold the result files')
    parser.add_argument('-a', '--adapters', metavar='/left_adapters.fasta',
                        required=True,
                        help='Fasta file with sequence of adapters to trim at the beginning of reads')
    parser.add_argument('-t', '--threads', metavar=str(max_cpu),
                        required=False,
                        type=int, default=max_cpu,
                        help='Number of CPU. Default is maximum CPU available({})'.format(max_cpu))
    parser.add_argument('-m', '--memory', metavar=str(max_cpu),
                        required=False,
                        type=int, default=max_mem,
                        help='Memory in GB. Default is 85% of total memory ({})'.format(max_mem))
    parser.add_argument('-pe', '--paired-end',
                        required=False,
                        action='store_true',  # implies default=False
                        help='Input fastq files are paired-end')
    parser.add_argument('-se', '--single-end',
                        required=False,
                        action='store_true',  # implies default=False
                        help='Input fastq files are single-end')

    # Get the arguments into an object
    arguments = parser.parse_args()

    GBS(arguments)
