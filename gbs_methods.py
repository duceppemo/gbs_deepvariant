#!/usr/local/env python3

from glob import glob
import subprocess
from multiprocessing import cpu_count
from psutil import virtual_memory
import sys
import gzip
import os
import pathlib
from concurrent import futures
import pysam
from natsort import natsort as ns
import pandas as pd
from deepvariant import DeepVariant


class SampleObject(object):
    def __init__(self, file_path):
        # Create seq object with its attributes
        self.file_path = file_path


class AdapterObject(object):
    def __init__(self, direction, seq):
        # Create seq object with its attributes
        self.direction = direction
        self.seq = seq


class Methods(object):
    accepted_extensions = ['.fq', '.fastq', '.fq.gz', '.fastq.gz']

    @staticmethod
    def get_files(in_folder, sample_dict):
        # Look for input sequence files recursively
        for root, directories, filenames in os.walk(in_folder):
            for filename in filenames:
                if filename.endswith(tuple(Methods.accepted_extensions)):  # accept a tuple or string
                    file_path = os.path.join(root, filename)
                    sample = filename.split('.')[0]
                    try:
                        if sample_dict[sample]:
                            sample_dict[sample].file_list.append(file_path)
                    except KeyError:
                        sample_dict[sample] = SampleObject([file_path])

    @staticmethod
    def check_folder_empty(folder):
        status = 0
        # List content of folder
        dir_content = os.listdir(folder)
        # if folder is not empty and all files have the accepted extensions
        test_file_ext = [x.endswith(tuple(Methods.accepted_extensions)) for x in dir_content]
        if dir_content and any(test_file_ext):
            status = 1
        return status

    @staticmethod
    def check_cpus(requested_cpu):
        total_cpu = cpu_count()

        if requested_cpu:
            if requested_cpu > total_cpu:
                requested_cpu = total_cpu
                sys.stderr.write("Number of threads was set higher than available CPUs ({})".format(total_cpu))
                sys.stderr.write("Number of threads was set to {}".format(requested_cpu))
        else:
            requested_cpu = total_cpu

        return requested_cpu

    @staticmethod
    def check_mem(requested_mem):
        max_mem = int(virtual_memory().total * 0.85 / 1000000000)  # in GB
        if requested_mem:
            if requested_mem > max_mem:
                requested_mem = max_mem
                sys.stderr.write("Requested memory was set higher than available system memory ({})".format(max_mem))
                sys.stderr.write("Memory was set to {}".format(requested_mem))
        else:
            requested_mem = max_mem

        return requested_mem

    @staticmethod
    def create_output_folders(output):
        """
        Create the output folder and subfolders
        :param output: string; absolute path to output folder
        :return:
        """
        # Create output folder is it does not exist
        pathlib.Path(output).mkdir(parents=True, exist_ok=True)

    @staticmethod
    def parse_adapters(adapter_file, adapter_dict):
        with open(adapter_file, 'r') as f:
            adapter_name = ''
            for line in f:
                line = line.rstrip()
                if not line:
                    continue
                if line.startswith('>'):
                    adapter_name = '-'.join(line.split('-')[:-1])[1:]  # drop the leading '>' and trailing '-R' or '-L'
                    if line.endswith('-L'):
                        adapter_dict[adapter_name] = AdapterObject('left', '')
                    elif line.endswith('-R'):
                        adapter_dict[adapter_name] = AdapterObject('right', '')
                    else:
                        raise Exception('Make sure your header has "-L" or "-R" at the end')
                elif adapter_name:
                    adapter_dict[adapter_name].seq = line
                else:
                    raise Exception('Something went wrong parsing the adapter file')

    @staticmethod
    def trim_both(fastq, output_folder, adapter_dict, mem, cpu):
        out_file = output_folder + os.path.basename(fastq)
        left_adapters = ','.join(y.seq for x, y in adapter_dict.items() if y.direction == 'left')[:-1]
        right_adapters = ','.join(y.seq for x, y in adapter_dict.items() if y.direction == 'right')[:-1]

        left_trim = ['bbduk.sh', '-Xmx{}g'.format(mem),
                     'threads={}'.format(cpu),
                     'in={}'.format(fastq),
                     'interleaved=f',
                     'editdistance=1',  # 1 mismatch allowed
                     'ktrim=l', 'k={}'.format(17), 'mink=15', 'rcomp=f',
                     'literal={}'.format(left_adapters),
                     'out=stdout.fq']
        right_trim = ['bbduk.sh', '-Xmx{}g'.format(mem),
                      'overwrite=t',
                      'threads={}'.format(cpu),
                      'in=stdin.fq',
                      'interleaved=f',
                      'editdistance=1',  # 1 mismatch allowed
                      'ktrim=r', 'k={}'.format(17), 'mink=15', 'rcomp=f',
                      'literal={}'.format(right_adapters),
                      'out={}'.format(out_file)]
        p1 = subprocess.Popen(left_trim, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        p2 = subprocess.Popen(right_trim, stdin=p1.stdout, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        p1.stdout.close()
        p2.communicate()

    @staticmethod
    def trim_left(fastq, output_folder, adapter_dict, mem, cpu, sample_dict):
        out_file = output_folder + os.path.basename(fastq)
        left_adapters = ','.join(y.seq for x, y in adapter_dict.items() if y.direction == 'left')

        left_trim = ['bbduk.sh', '-Xmx{}g'.format(mem),
                     'threads={}'.format(cpu),
                     'in={}'.format(fastq),
                     'editdistance=1',  # 1 mismatch allowed
                     'ktrim=l', 'k={}'.format(17), 'mink=15', 'rcomp=f',
                     'literal={}'.format(left_adapters),
                     'minlength=64',
                     'out={}'.format(out_file)]

        subprocess.run(left_trim, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

    @staticmethod
    def trim_right(fastq, output_folder, adapter_dict, mem, cpu, sample_dict):
        out_file = output_folder + os.path.basename(fastq)
        right_adapters = ','.join(y.seq for x, y in adapter_dict.items() if y.direction == 'right')

        right_trim = ['bbduk.sh', '-Xmx{}g'.format(mem),
                      'threads={}'.format(cpu),
                      'in={}'.format(fastq),
                      'editdistance=1',  # 1 mismatch allowed
                      'ktrim=r', 'k={}'.format(17), 'mink=15', 'rcomp=f',
                      'literal={}'.format(right_adapters),
                      'minlength=64',
                      'out={}'.format(out_file)]

        subprocess.run(right_trim, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

    @staticmethod
    def parallel_trim_reads(trim_function, adapter_dict, sample_dict, output_folder, mem, cpu):
        print('Trimming reads...')
        Methods.create_output_folders(output_folder)

        with futures.ThreadPoolExecutor(max_workers=cpu/8) as executor:
            # fastq, output_folder, adapter_dict, mem, cpu
            args = ((info.file_path[0], output_folder, adapter_dict, mem, cpu)
                    for sample, info in sample_dict.items())
            for results in executor.map(lambda x: trim_function(*x), args):
                pass

    @staticmethod
    def index_bowtie2(ref, prefix, cpu):
        print('Indexing reference genome...')
        cmd = ['bowtie2-build', '--threads', str(cpu), ref, prefix]
        subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

    @staticmethod
    def map_bowtie2(ref_index, fastq_file, cpu, output_bam):
        sample = os.path.basename(fastq_file).replace('.fastq.gz', '')
        rg = '@RG\tID:{}\tSM:{}'.format(sample, sample)
        bowtie2_align_cmd = ['bowtie2',
                             '--very-sensitive-local',
                             '--xeq ',
                             '-x', ref_index,
                             '-U', fastq_file,
                             '--threads', str(cpu),
                             '--rg', rg]
        samtools_view_cmd = ['samtools', 'view',
                             '-@', str(cpu),
                             '-F', '4',
                             '-b', '-h', '-']
        samtools_sort_cmd = ['samtools', 'sort',
                             '-@', str(cpu),
                             '-o', output_bam,
                             '-']
        samtools_index_cmd = ['samtools', 'index',
                              output_bam]

        p1 = subprocess.Popen(bowtie2_align_cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        p2 = subprocess.Popen(samtools_view_cmd, stdin=p1.stdout, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        p1.stdout.close()
        p3 = subprocess.Popen(samtools_sort_cmd, stdin=p2.stdout, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        p2.stdout.close()
        p3.communicate()

        # index bam file
        subprocess.run(samtools_index_cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

    @staticmethod
    def parallel_map_bowtie2(output_folder, ref, sample_dict, cpu):
        print('Mapping reads...')
        Methods.create_output_folders(output_folder)

        with futures.ThreadPoolExecutor(max_workers=cpu/8) as executor:
            args = ((ref, info.file_path[0], 8, output_folder + '/' + sample + '.bam')
                    for sample, info in sample_dict.items())
            for results in executor.map(lambda x: Methods.map_bowtie2(*x), args):
                pass

    @staticmethod
    def index_samtools(fasta):
        pysam.faidx(fasta)

    @staticmethod
    def list_to_file(my_list, output_file):
        with open(output_file, 'wt') as f:
            for l in my_list:
                f.write('{}\n'.format(l))

    @staticmethod
    def call_variants(ref, bam_list_file, vcf, cpu):
        pileup_cmd = ['bcftools', 'mpileup',
                      '-b', bam_list_file,
                      '-O', 'u',
                      '-f', ref,
                      '--threads', str(cpu)]
        call_cmd = ['bcftools', 'call',
                    '-Ov',
                    '-m', '-v',
                    '-o', vcf,
                    '--threads', str(cpu)]
        p1 = subprocess.Popen(pileup_cmd, stdout=subprocess.PIPE)  # , stderr=subprocess.DEVNULL)
        p2 = subprocess.Popen(call_cmd, stdin=p1.stdout, stdout=subprocess.PIPE)  # , stderr=subprocess.DEVNULL)
        p1.stdout.close()
        p2.communicate()

    @staticmethod
    def activate_deepvariant():
        pass

    @staticmethod
    def call_variants_deepvariant(container, ref, bam_list, vcf, cpu, output_dir):
        # Assume that deepvariant docker image installed and user is part of docker group so don't need to run sud
        # Check if container is available
        container_name = 'deepvariant_gpu'
        if not container:  # TODO -> replace with actual command to see if container is available or not
            DeepVariant.get_container(container_name,
                                      'gcr.io/deepvariant-docker/{}:{}'.format(
                                          container_name, DeepVariant.BIN_VERSION), output_dir)
        # Check if container is running
        elif not container:  #TODO -> container is present, but not started
            DeepVariant.start_container(container_name)

        # Create folder structure for output files
        DeepVariant.set_folders(output_dir)

        #Download model files
        DeepVariant.download_model(output_dir)

    @staticmethod
    def fix_vcf(vcf_in, vcf_out):
        with open(vcf_out, 'w') as out_fh:
            with open(vcf_in, 'r') as in_fh:
                for line in in_fh:
                    line = line.rstrip()
                    if not line:
                        continue
                    if line.startswith('##contig'):
                        continue
                    elif line.startswith('##'):
                        out_fh.write('{}\n'.format(line))
                    elif line.startswith('#CHROM'):
                        field_list = line.split('\t')
                        sample_list = field_list[9:]
                        sample_list = ['.'.join(os.path.basename(file_name).split('.')[:-1])
                                       for file_name in sample_list]
                        out_fh.write('{}\t{}\n'.format('\t'.join(field_list[:9]), '\t'.join(sample_list)))
                    else:
                        field_list = line.split('\t')
                        # add an ID
                        field_list[2] = '{}:{}'.format(field_list[0], field_list[1])
                        out_fh.write('{}\n'.format('\t'.join(field_list)))

    @staticmethod
    def sort_vcf(vcf_in, vcf_out):
        # https://gist.github.com/slowkow/6215557
        lines_to_skip = 0
        headers = list()
        with open(vcf_out, 'w') as out_fh:
            with open(vcf_in, 'r') as in_fh:
                for line in in_fh:
                    line = line.rstrip()
                    if not line:
                        continue
                    if line.startswith('#'):
                        lines_to_skip += 1
                        if line.startswith('#CHROM'):
                            headers = line.split('\t')
                        else:
                            out_fh.write('{}\n'.format(line))
                    else:
                        break
            # Parse VCF sample info to dataframe
            df = pd.read_table(vcf_in, sep='\t', skiprows=lines_to_skip, names=headers)
            # Sort dataframe
            df.sort_values(by=['#CHROM', 'POS'], ascending=True)
            df.to_csv(out_fh, sep='\t', header=True, index=False)

    @staticmethod
    def filter_variants(vcf_in, vcf_out):
        pass

    @staticmethod
    def impute_missing_snps():
        pass

    @staticmethod
    def LD_filter(vcf_in, vcf_out):
        pass

    @staticmethod
    def gwas(vcf_in, phenotypes, output_file):
        pass

    @staticmethod
    def make_manhattan_plot(input_file, output_file):
        pass

