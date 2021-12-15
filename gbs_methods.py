#!/usr/local/env python3

# docker pull

# mamba create -n gbs_deepvariant
# mamba activate gbs_deepvariant
# mamba install -c bioconda bbmap samtools bowtie2 pysam bcftools vcftools -y
# mamba install -c conda-forge psutil -y


import random
import subprocess
from multiprocessing import cpu_count
from psutil import virtual_memory
import sys
import os
import pathlib
from concurrent import futures
import pysam
import pandas as pd
import io
import gzip

import gbs_methods


class VcfObject(object):
    def __init__(self, chrom, pos, ident, ref, alt, qual, filt, info, form, geno):
        # Create object with its attributes
        self.chrom = chrom
        self.pos = pos
        self.ident = ident
        self.ref = ref
        self.alt = alt
        self.qual = qual
        self.filt = filt
        self.info = info
        self.form = form
        self.geno = geno


class AdapterObject(object):
    def __init__(self, direction, seq):
        # Create seq object with its attributes
        self.direction = direction
        self.seq = seq


class Methods(object):
    accepted_extensions = ['.fq', '.fastq', '.fq.gz', '.fastq.gz']

    vcf_header = ['##fileformat=VCFv4.2',
                  '##FILTER=<ID=PASS,Description="All filters passed">'
                  '##FILTER=<ID=RefCall,Description="Genotyping model thinks this site is reference.">',
                  '##FILTER=<ID=LowQual,Description="Confidence in this variant being real is below calling threshold.">',
                  '##FILTER=<ID=NoCall,Description="Site has depth=0 resulting in no call.">',
                  '##INFO=<ID=END,Number=1,Type=Integer,Description="End position (for use with symbolic alleles)">',
                  '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                  '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Conditional genotype quality">',
                  '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">',
                  '##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="Minimum DP observed within the GVCF block.">',
                  '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Read depth for each allele">',
                  '##FORMAT=<ID=VAF,Number=A,Type=Float,Description="Variant allele fractions.">',
                  '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled genotype likelihoods rounded to the closest integer">',
                  '##FORMAT=<ID=MED_DP,Number=1,Type=Integer,Description="Median DP observed within the GVCF block rounded to the nearest integer.">',
                  '##DeepVariant_version=1.2.0',
                  '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT']

    @staticmethod
    def make_folder(folder):
        """
        Create output folder.
        :param folder: string. Output folder path.
        :return:
        """
        # Will create parent directories if don't exist and will not return error if already exists
        pathlib.Path(folder).mkdir(parents=True, exist_ok=True)

    @staticmethod
    def get_files(in_folder, sample_dict):
        # Look for input sequence files recursively
        for root, directories, filenames in os.walk(in_folder):
            for filename in filenames:
                if filename.endswith(tuple(Methods.accepted_extensions)):  # accept a tuple or string
                    file_path = os.path.join(root, filename)
                    sample = filename.split('.')[0].split('_R1')[0].split('_R2')[0]
                    if '_R1' in filename:
                        sample_dict[sample].insert(0, file_path)
                    elif '_R2' in filename:
                        sample_dict[sample].insert(1, file_path)

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
    def create_folders(output):
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
    def parallel_trim_reads(trim_function, adapter_dict, sample_dict, output_folder, mem, cpu, p):
        print('Trimming reads...')
        Methods.create_folders(output_folder)

        with futures.ThreadPoolExecutor(max_workers=p) as executor:
            # fastq, output_folder, adapter_dict, mem, cpu
            args = ((info.file_path[0], output_folder, adapter_dict, mem, int(cpu/p))
                    for sample, info in sample_dict.items())
            for results in executor.map(lambda x: trim_function(*x), args):
                pass

    @staticmethod
    def index_bowtie2(ref, prefix, cpu):
        print('Indexing reference genome...')
        cmd = ['bowtie2-build', '--threads', str(cpu), ref, prefix]
        subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

    @staticmethod
    def map_bowtie2_se(ref, fastq_file, cpu, output_bam):
        sample = os.path.basename(fastq_file).replace('.fastq.gz', '')
        ref_index = '.'.join(ref.split('.')[:-1])

        print('\t{}'.format(sample))
        rg = '@RG\tID:{}\tSM:{}'.format(sample, sample)
        bowtie2_align_cmd = ['bowtie2',
                             # '--very-sensitive-local',
                             '--xeq ',
                             '-x', ref_index,
                             '-U', fastq_file,
                             '--threads', str(cpu),
                             '--rg', rg]
        samtools_view_cmd = ['samtools', 'view',
                             '-@', str(cpu),
                             '-F', '4',
                             '-C', '-h',
                             '-T', ref,
                             '-']
        samtools_sort_cmd = ['samtools', 'sort',
                             '-@', str(cpu),
                             '-o', output_bam,
                             '-']
        # samtools can only index chromosomes up to 512M bp.
        samtools_index_cmd = ['samtools', 'index',
                              '-c',
                              output_bam]

        p1 = subprocess.Popen(bowtie2_align_cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        p2 = subprocess.Popen(samtools_view_cmd, stdin=p1.stdout, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        p1.stdout.close()
        p3 = subprocess.Popen(samtools_sort_cmd, stdin=p2.stdout, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        p2.stdout.close()
        p3.communicate()

        # index cram file
        subprocess.run(samtools_index_cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

    @staticmethod
    def parallel_map_bowtie2_se(output_folder, ref, sample_dict, cpu, p):
        # Index reference genome if not already done
        ref_index = '.'.join(ref.split('.')[:-1])
        if not os.path.exists(ref_index + '.1.bt2') and not os.path.exists(ref_index + '.1.bt2l'):
            Methods.index_bowtie2(ref, ref_index, cpu)

        print('Mapping reads...')
        Methods.create_folders(output_folder)

        with futures.ThreadPoolExecutor(max_workers=p) as executor:
            args = ((ref, path_list[0], int(cpu/p), output_folder + '/' + sample + '.cram')
                    for sample, path_list in sample_dict.items())
            for results in executor.map(lambda x: Methods.map_bowtie2_se(*x), args):
                pass

    @staticmethod
    def map_bowtie2_pe(ref, r1, r2, cpu, output_bam):
        sample = os.path.basename(r1).replace('.fastq.gz', '').split('_R1')[0]  # TODO: search and replace allowed_ext
        ref_index = '.'.join(ref.split('.')[:-1])

        print('\t{}'.format(sample))
        # rg = '@RG\tID:{}\tSM:{}'.format(sample, sample)
        bowtie2_align_cmd = ['bowtie2',
                             '--very-sensitive-local',
                             '--xeq ',
                             '-x', ref_index,
                             '-1', r1,
                             '-2', r2,
                             '--threads', str(cpu),
                             '--rg-id', '{}.{}'.format(sample, random.randrange(0, 1000000, 1)),
                             '--rg', 'ID:{}'.format(sample),
                             '--rg', 'SM:{}'.format(sample)]
        samtools_view_cmd = ['samtools', 'view',
                             '-@', str(cpu),
                             '-F', '4',
                             '-C', '-h',
                             '-T', ref,
                             '-']
        samtools_sort_cmd = ['samtools', 'sort',
                             '-@', str(cpu),
                             '-o', output_bam,
                             '-']
        # samtools can only index chromosomes up to 512M bp.
        samtools_index_cmd = ['samtools', 'index',
                              '-c',
                              output_bam]

        p1 = subprocess.Popen(bowtie2_align_cmd, stdout=subprocess.PIPE)#, stderr=subprocess.DEVNULL)
        p2 = subprocess.Popen(samtools_view_cmd, stdin=p1.stdout, stdout=subprocess.PIPE)#, stderr=subprocess.DEVNULL)
        p1.stdout.close()
        p3 = subprocess.Popen(samtools_sort_cmd, stdin=p2.stdout, stdout=subprocess.PIPE)#, stderr=subprocess.DEVNULL)
        p2.stdout.close()
        p3.communicate()

        # index cram file
        subprocess.run(samtools_index_cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

    @staticmethod
    def parallel_map_bowtie2_pe(output_folder, ref, sample_dict, cpu, p):
        # Index reference genome if not already done
        ref_index = '.'.join(ref.split('.')[:-1])
        if not os.path.exists(ref_index + '.1.bt2') and not os.path.exists(ref_index + '.1.bt2l'):
            Methods.index_bowtie2(ref, ref_index, cpu)

        print('Mapping reads...')
        Methods.create_folders(output_folder)

        with futures.ThreadPoolExecutor(max_workers=p) as executor:
            args = ((ref, path_list[0], path_list[1], int(cpu/p), output_folder + '/' + sample + '.cram')
                    for sample, path_list in sample_dict.items())
            for results in executor.map(lambda x: Methods.map_bowtie2_pe(*x), args):
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
    def move_files(input_folder, output_folder, ext):
        to_move_list = list()
        for root, directories, filenames in os.walk(input_folder):
            for filename in filenames:
                if filename.endswith(ext):  # accept a tuple or string
                    source_file = os.path.join(root, filename)
                    destination_file = os.path.join(output_folder, filename)
                    if source_file.endswith(ext):
                        os.replace(source_file, destination_file)

    @staticmethod
    def vcf_to_pandas(vcf):
        # https://gist.github.com/dceoy/99d976a2c01e7f0ba1c813778f9db744
        with open(vcf, 'r') as f:
            lines = [l for l in f if not l.startswith('##')]
        return pd.read_csv(io.StringIO(''.join(lines)), dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str,
                                                               'ALT': str, 'QUAL': str, 'FILTER': str, 'INFO': str},
                           sep='\t').rename(columns={'#CHROM': 'CHROM'})

    @staticmethod
    def vcf_to_dict(vcf_file):
        sample = '.'.join(os.path.basename(vcf_file).split('.')[:-1])
        vcf_dict = dict()
        with gzip.open(vcf_file, 'rb') if vcf_file.endswith('gz') else open(vcf_file, 'r') as f:
            for line in f:
                if vcf_file.endswith('gz'):
                    line = line.decode()
                line = line.rstrip()
                if line.startswith('#'):
                    continue
                else:
                    (chrom, pos, ident, ref, alt, qual, filt, info, form, geno) = line.split('\t')
                    k = '{}:{}'.format(chrom, pos)
                    try:
                        vcf_dict[k]
                    except KeyError:
                        vcf_dict[k] = dict()

                    vcf_dict[k]['REF'] = ref
                    vcf_dict[k]['ALT'] = alt
                    vcf_dict[k]['QUAL'] = qual
                    vcf_dict[k]['FILTER'] = filt
                    vcf_dict[k]['INFO'] = info
                    vcf_dict[k]['FORMAT'] = form
                    vcf_dict[k]['GENO'] = [(sample, geno)]

        return vcf_dict

    @staticmethod
    def merge_vcf(input_folder, output_vcf):
        vcf_list = list()
        for root, directories, filenames in os.walk(input_folder):
            for filename in filenames:
                if filename.endswith('.vcf'):  # accept a tuple or string
                    vcf_list.append(os.path.join(root, filename))

        sample_list = ['.'.join(os.path.basename(x).split('.')[:-1]) for x in vcf_list]
        # main_df = pd.DataFrame(columns=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'])

        # Parse VCF files
        main_dict = dict()
        for vcf in vcf_list:
            # Convert VCF file to dictionary
            vcf_dict = Methods.vcf_to_dict(vcf)
            for k, info in vcf_dict.items():
                try:
                    main_dict[k]
                except KeyError:
                    main_dict[k] = dict()
                if not main_dict[k]:
                    main_dict[k] = info
                else:
                    main_dict[k]['GENO'] += info['GENO']
            # sample_df = Methods.vcf_to_pandas(vcf)
            # main_df = main_df.append(sample_df, sort=False)

        # Get header
        my_vcf = vcf_list[0]
        header_list = list()
        with gzip.open(my_vcf, 'rb') if my_vcf.endswith('gz') else open(my_vcf, 'r') as f:
            for line in f:
                if my_vcf.endswith('gz'):
                    line = line.decode()
                line = line.rstrip()
                if line.startswith('##'):
                    header_list.append(line)
        header_list.append('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT')

        # Ouput merged VCF files
        with open(output_vcf, 'w') as f:
            # Write VCF header
            f.write('{}\t{}\n'.format('\n'.join(header_list), '\t'.join(sample_list)))

            # Merge the genotype information for all the samples
            for k, info in main_dict.items():
                geno = ['./.'] * len(sample_list)  # Assign missing value by default
                for geno_tuple in info['GENO']:
                    i = sample_list.index(geno_tuple[0])
                    geno[i] = geno_tuple[1]

                chrom, pos = k.split(':')
                ident = '.'
                ref = info['REF']
                alt = info['ALT']
                qual = info['QUAL']
                filt = info['FILTER']
                inf = info['INFO']
                form = info['FORMAT']
                f.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t\n'.format(chrom, pos, ident, ref, alt, qual, filt,
                                                                            inf, form, '\t'.join(geno)))

        # Sort merged VCF file
        print('Sorting merged VCF file...')
        Methods.sort_vcf(output_vcf, output_vcf + '.sorted')
        os.replace(output_vcf + '.sorted', output_vcf)

        # Sort df by 1) CHROM and by 2) POS
        # main_df = main_df.sort_values(['CHROM', 'POS'], ascending=[True, True])
        # Combine rows if value of CHROM and POS is the same
        # Write dataframe to file
        # main_df.to_csv(output_vcf, sep='\t', index=False, na_rep='./.')

    @staticmethod
    def sort_vcf(vcf_in, vcf_out):
        cmd = ['bcftools', 'sort',
               '-o', vcf_out,
               vcf_in]
        subprocess.run(cmd)

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
                        sample_list = ['.'.join(
                            os.path.basename(file_name).split('.')[:-1]) if '.gz' not in file_name else '.'.join(
                            os.path.basename(file_name).split('.')[:-2]) for file_name in sample_list]
                        out_fh.write('{}\t{}\n'.format('\t'.join(field_list[:9]), '\t'.join(sample_list)))
                    else:
                        field_list = line.split('\t')
                        # add an ID
                        field_list[2] = '{}:{}'.format(field_list[0], field_list[1])
                        out_fh.write('{}\n'.format('\t'.join(field_list)))

    # @staticmethod
    # def sort_vcf(vcf_in, vcf_out):
    #     # https://gist.github.com/slowkow/6215557
    #     lines_to_skip = 0
    #     headers = list()
    #     with open(vcf_out, 'w') as out_fh:
    #         with open(vcf_in, 'r') as in_fh:
    #             for line in in_fh:
    #                 line = line.rstrip()
    #                 if not line:
    #                     continue
    #                 if line.startswith('#'):
    #                     lines_to_skip += 1
    #                     if line.startswith('#CHROM'):
    #                         headers = line.split('\t')
    #                     else:
    #                         out_fh.write('{}\n'.format(line))
    #                 else:
    #                     break
    #         # Parse VCF sample info to dataframe
    #         df = pd.read_table(vcf_in, sep='\t', skiprows=lines_to_skip, names=headers)
    #         # Sort dataframe
    #         df.sort_values(by=['#CHROM', 'POS'], ascending=True)
    #         df.to_csv(out_fh, sep='\t', header=True, index=False, na_rep='./.')

    @staticmethod
    def filter_variants(vcf_in, vcf_out):
        # https://vcftools.github.io/man_latest.html
        cmd = ['vcftools',
               '--vcf', vcf_in,
               '--out', vcf_out,
               '--remove-indels',
               '--maf', str(0.05),
               '--min-alleles', str(2),
               '--max-alleles', str(2),
               '--minDP', str(5),  # vs '--min-meanDP'
               '--minQ', str(20),
               '--remove-filtered-geno-all',
               '--recode']

        subprocess.run(cmd)

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

