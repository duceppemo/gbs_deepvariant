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


from argparse import ArgumentParser
from multiprocessing import cpu_count
from psutil import virtual_memory
import os
from glob import glob
from gbs_methods import Methods
from deepvariant import DeepVariant
from collections import defaultdict


class GBS(object):
    def __init__(self, args):
        # Command line arguments
        self.input = args.input
        self.out_folder = args.output
        self.cpu = args.threads
        self.parallel = args.parallel
        self.mem = args.memory
        self.ref = args.reference
        if args.paired_end:
            self.read_type = 'pe'
        elif args.single_end:
            self.read_type = 'se'
        elif (args.paired_end and args.single_end)\
                or (not args.paired_end and not args.single_end):
            raise Exception('Please choose between "-se" and "-pe"\n')
        self.ion = args.ion_torrent

        # Read length auto removal
        self.size = args.size

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

        # Step completion report files
        done_trimming = self.out_folder + '/done_trimming'
        done_mapping = self.out_folder + '/done_mapping'
        done_deepvariant = self.out_folder + '/done_deepvariant'
        done_merging = self.out_folder + '/done_merging'
        done_converting = self.out_folder + '/done_converting'

        # Output folders to create
        trimmed_folder = self.out_folder + '/1_trimmed/'
        mapped_folder = self.out_folder + '/2_mapped/'
        deepvariant_folder = self.out_folder + '/3_deepvariant/'
        gvcf_folder = self.out_folder + '/tmp_gvcf/'
        merged_gvcf_folder = self.out_folder + '/4_merged_gvcf/'
        merged_vcf_folder = self.out_folder + '/5_merged_vcf/'
        filtered_folder = self.out_folder + '/6_filtered/'

        # Get input files and place info in dictionary
        self.sample_dict['raw'] = Methods.get_files(self.input)

        # Trim reads
        if not os.path.exists(done_trimming):
            # Auto size selection
            if self.size == 'auto':
                # Find the best read length for trimming
                print('Auto detecting min read size to keep...')
                df = Methods.parallel_read_length_dist(self.sample_dict['raw'], self.cpu)

                # Create plot
                Methods.make_folder(read_length)
                fig1 = px.line(df, x=df.index, y='Count')
                with open(read_length + '/' + 'first_figure.html', 'w') as fh:
                    fh.write(fig1.to_html(full_html=False, include_plotlyjs='cdn'))

                # Detect peak
                self.size = Methods.find_peak(df, read_length)

            # Trim reads
            if self.ion:  # if Ion Torrent reads
                print('Processing IonTorrent reads...')
                Methods.parallel_trim_reads(Methods.trim_iontorrent, self.sample_dict['raw'], trimmed_folder,
                                            self.cpu, self.parallel, self.size)
                Methods.flag_done(done_trimming)
            else:  # if Illumina reads
                if self.read_type == 'se':  # if single-end
                    print('Processing Illumina single-end reads...')
                    Methods.parallel_trim_reads(Methods.trim_illumina_se, self.sample_dict['raw'],
                                                trimmed_folder, self.cpu, self.parallel, self.size)
                    Methods.flag_done(done_trimming)
                else:  # elif self.read_type == 'pe':  # if paired-end
                    print('Processing Illumina paired-end reads...')
                    Methods.parallel_trim_reads(Methods.trim_illumina_pe, self.sample_dict['raw'],
                                                trimmed_folder, self.cpu, self.parallel, self.size)
                    Methods.flag_done(done_trimming)
                os.remove('fastp.json')  # fastp artefact
        else:
            print('Skipping trimming. Already done.')

        # Update sample_dict after trimming
        self.sample_dict['trimmed'] = Methods.get_files(trimmed_folder)

        # if self.adapters:
        #     Methods.parse_adapters(self.adapters, adapter_dict)
        #     if any(info.direction == 'left' for ident, info in adapter_dict.items())\
        #             and not any(info.direction == 'right' for ident, info in adapter_dict.items()):
        #         Methods.parallel_trim_reads(Methods.trim_left, adapter_dict, self.sample_dict,
        #                                     self.out_folder + '/trimmed/', self.mem, self.cpu, self.parallel)
        #     elif not any(info.direction == 'left' for ident, info in adapter_dict.items())\
        #             and any(info.direction == 'right' for ident, info in adapter_dict.items()):
        #         Methods.parallel_trim_reads(Methods.trim_right, adapter_dict, self.sample_dict,
        #                                     self.out_folder + '/trimmed/', self.mem, self.cpu, self.parallel)
        #     elif any(info.direction == 'left' for ident, info in adapter_dict.items())\
        #             and any(info.direction == 'right' for ident, info in adapter_dict.items()):
        #         Methods.parallel_trim_reads(Methods.trim_both, adapter_dict, self.sample_dict,
        #                                     self.out_folder + '/trimmed/', self.mem, self.cpu, self.parallel)

        # Map reads
        if not os.path.exists(done_mapping):
            if self.read_type == 'se':
                Methods.parallel_map_bowtie2_se(mapped_folder, self.ref, self.sample_dict['trimmed'],
                                                self.cpu, self.parallel)
            elif self.read_type == 'pe':
                Methods.parallel_map_bowtie2_pe(mapped_folder, self.ref, self.sample_dict['trimmed'],
                                                self.cpu, self.parallel)
            Methods.flag_done(done_mapping)
        else:
            print('Skipping mapping. Already done.')

        # Call variants
        if not os.path.exists(done_deepvariant):
            if not os.path.exists(self.ref + '.fai'):
                Methods.index_samtools(self.ref)

            # List bam files into a text file, one per line
            cram_list = glob(mapped_folder + '*.cram')
            Methods.list_to_file(cram_list, self.out_folder + '/cram.list')

            print('Calling vairants with DeepVariant v1.2.0...')
            DeepVariant.call_variants(self.ref, cram_list, deepvariant_folder, self.cpu)

            # Move gVCF files
            Methods.create_folders(gvcf_folder)
            Methods.move_files(deepvariant_folder, gvcf_folder, '.gvcf')

            # BGzip gvcf files
            Methods.parallel_bgzip_files(gvcf_folder, self.cpu, self.parallel)

            # Done flag
            Methods.flag_done(done_deepvariant)
        else:
            print('Skipping variant calling. Already done.')

        # Merge gVCF files
        merged_gvcf = merged_gvcf_folder + 'merged.gvcf.gz'
        if not os.path.exists(done_merging):
            print('Merging gVCF files...')
            Methods.create_folders(merged_gvcf_folder)
            Methods.merge_gvcf(self.ref, gvcf_folder, merged_gvcf, self.cpu)  # Merge with bcftools
            Methods.parallel_bgzip_files(merged_gvcf_folder, self.cpu, self.parallel)  # BGzip gvcf files
            Methods.flag_done(done_merging)
        else:
            print('Skipping gVCF merging. Already done.')

        # Convert gVCF to VCF
        merged_vcf = merged_vcf_folder + 'merged.vcf'
        fixed_vcf = merged_vcf_folder + 'merged_fixed.vcf'
        if not os.path.exists(done_converting):
            print('Converting gVCF to VCF...')
            Methods.create_folders(merged_vcf_folder)
            # Methods.convert_gvcf_to_vcf(self.ref, merged_gvcf, merged_vcf, self.cpu)
            Methods.fix_merged_vcf(merged_vcf, fixed_vcf)
            Methods.flag_done(done_converting)
        else:
            print('Skipping gVCF conversion to VCF. Already done.')

        # Filtering VCF file
        filtered_vcf = filtered_folder + 'filtered'
        Methods.create_folders(filtered_folder)
        Methods.filter_variants(fixed_vcf, filtered_vcf)
        os.replace(filtered_folder + '/filtered.recode.vcf', filtered_folder + '/filtered.vcf')


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
                        required=False,
                        help='Fasta file with sequence of adapters to trim at the beginning of reads')
    parser.add_argument('-s', '--size', metavar='64|auto',
                        required=False, default=64,
                        help='Minimum read size to keep after trimming. '
                             'An experimetial auto size selection (use "auto" as argument) is also available. '
                             'It is based on highest peak detection after plotting read length distribution. '
                             'Experimental. Default is 64. Optional.')
    parser.add_argument('-t', '--threads', metavar=str(max_cpu),
                        required=False,
                        type=int, default=max_cpu,
                        help='Number of CPU. Default is maximum CPU available({})'.format(max_cpu))
    parser.add_argument('-p', '--parallel', metavar='4',
                        required=False,
                        type=int, default=4,
                        help='Number of samples to run in parallel')
    parser.add_argument('-m', '--memory', metavar=str(max_mem),
                        required=False,
                        type=int, default=max_mem,
                        help='Memory in GB. Default is 85%% of total memory ({})'.format(max_mem))
    parser.add_argument('-ion', '--ion-torrent',
                        required=False,
                        action='store_true',  # implies default=False
                        help='Reads are from IonTorrent (different trim parameters). Default is Illumina.')
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
