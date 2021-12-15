#!/usr/local/env python3

# from gbs_methods import Methods
import subprocess
import os
from pathlib import Path


class DeepVariant(object):
    @staticmethod
    def call_variants(ref, cram_list, deepvariant_folder, cpu):
        ref_folder = os.path.dirname(ref)
        ref_name = os.path.basename(ref)
        cram_folder = os.path.dirname(cram_list[0])

        for cram in cram_list:
            cram_name = '.'.join(os.path.basename(cram).split('.')[:-1])  # only remove what is after the last "."
            output_folder = deepvariant_folder + cram_name
            DeepVariant.run_deepvariant(ref_folder, ref_name, cram_folder, cram_name, output_folder, cpu)

    @staticmethod
    def run_deepvariant(ref_folder, ref_name, cram_folder, cram_name, output_folder, cpu):
        # Create output folder
        Path(output_folder).mkdir(parents=True, exist_ok=True)

        cmd = ['docker', 'run',
               '-v', '{}:/ref'.format(ref_folder),
               '-v', '{}:/input'.format(cram_folder),
               '-v', '{}:/output'.format(output_folder),
               'google/deepvariant:1.2.0',
               '/opt/deepvariant/bin/run_deepvariant',
               '--model_type={}'.format('WGS'),  # [WGS,WES,PACBIO,HYBRID_PACBIO_ILLUMINA]
               '--ref=/ref/{}'.format(ref_name),
               '--reads=/input/{}'.format(cram_name + '.cram'),
               '--sample_name={}'.format(cram_name),
               '--output_vcf=/output/{}.vcf'.format(cram_name),
               '--output_gvcf=/output/{}.gvcf'.format(cram_name),
               # '--call_variants_extra_args="use_openvino=true"',  # Optional for Intel CPUs reduces call_variants time
               '--num_shards={}'.format(cpu),  # use all cores for make_examples
               '--logging_dir=/output/logs',
               '--dry_run=false']

        subprocess.run(cmd)

    # def run_deepvariant(ref_folder, ref_name, cram_folder, cram_name, output_folder, cpu):
    #     # https://github.com/google/deepvariant/blob/r1.2/docs/deepvariant-quick-start.md#notes-on-gpu-image
    #     # https://github.com/google/deepvariant/blob/r1.2/scripts/install_nvidia_docker.sh
    #     # docker pull google/deepvariant:1.2.0-gpu
    #
    #     # Create output folder
    #     Path(output_folder).mkdir(parents=True, exist_ok=True)
    #
    #     cmd = ['docker', 'run',
    #            '-v', '{}:/ref'.format(ref_folder),
    #            '-v', '{}:/input'.format(cram_folder),
    #            '-v', '{}:/output'.format(output_folder),
    #            'google/deepvariant:1.2.0',
    #            '/opt/deepvariant/bin/run_deepvariant',
    #            '--model_type={}'.format('WGS'),  # [WGS,WES,PACBIO,HYBRID_PACBIO_ILLUMINA]
    #            '--ref=/ref/{}'.format(ref_name),
    #            '--reads=/input/{}'.format(cram_name + '.cram'),
    #            '--sample_name={}'.format(cram_name),
    #            '--output_vcf=/output/{}.vcf'.format(cram_name),
    #            '--output_gvcf=/output/{}.gvcf'.format(cram_name),
    #            # '--call_variants_extra_args="use_openvino=true"',  # Optional for Intel CPUs reduces call_variants time
    #            '--num_shards={}'.format(cpu),  # use all cores for make_examples
    #            '--logging_dir=/output/logs',
    #            '--dry_run=false']
    #
    #     subprocess.run(cmd)

