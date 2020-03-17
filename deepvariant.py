#!/usr/local/env python3

# from gbs_methods import Methods
import subprocess
import os


class DeepVariant(object):
    BIN_VERSION = "0.8.0"
    MODEL_VERSION = "0.8.0"

    MODEL_NAME = 'DeepVariant-inception_v3-{}+data-wgs_standard'.format(MODEL_VERSION)
    MODEL_HTTP_DIR = 'https://storage.googleapis.com/deepvariant/models/DeepVariant/{}/{}'.format(
        MODEL_VERSION, MODEL_NAME)
    DATA_HTTP_DIR = 'https://storage.googleapis.com/deepvariant/quickstart-testdata'

    @staticmethod
    def start_container(name):
        # Methods.create_output_folders(name)
        pass

    @staticmethod
    def get_container(name, address, DATA_DIR):
        cmd = ['docker', 'run',
               '--runtime=nvidia',
               '--name', name,
               '-v', '{}:/home/workspace'.format(DATA_DIR),
               address]
        subprocess.run(cmd)

    @staticmethod
    def set_folders(output_dir):
        for f in ['models', 'output', 'bam', 'log', 'ref']:
            DeepVariant.set_folders(output_dir + '/' + f)

    @staticmethod
    def download_model(output_folder):
        if not os.path.exists():
            pass
