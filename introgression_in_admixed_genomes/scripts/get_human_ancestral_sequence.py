#!/usr/bin/env python3
import argparse
import sys
from ftplib import FTP
import logging
import hashlib
import os


class EPO:
    def __init__(self, ftp_dir, epo_dir, force_overwrite=False):
        """
        Initialize EPO class
        :param ftp_dir: str, path to FTP directory from where to download EPO alignment files
        :param epo_dir:
        :param force_overwrite:
        """
        self.ftp_dir = ftp_dir
        if epo_dir:
            self.output_dir = epo_dir
        else:
            if ftp_dir.endswith('/'):
                self.output_dir = ftp_dir.split('//')[1].split('/')[-2].replace('.', '_')
            else:
                self.output_dir = ftp_dir.split('//')[1].split('/')[-1].replace('.', '_')
        self.overwrite = force_overwrite
        self.host = ftp_dir.split('//')[1].split('/')[0]
        logging.info(f'Connecting to host: {self.host}')
        self.ftp = FTP(self.host)
        self.ftp.login()
        self.directory = '/'.join(ftp_dir.split('//')[1].split('/')[1:])
        self.ftp.cwd(self.directory)
        logging.info(f'Changing working directory to: {self.directory}')
        self.files = self.ftp.nlst()
        if not os.path.isdir(self.output_dir):
            logging.info(f'Creating: {self.output_dir}')
            os.makedirs(self.output_dir)
        self.md5sums = self.get_md5sums

    def get_md5sums(self):
        """
        Download MD5SUM file from FTP directory and parse file to dict
        :return: dict, File with md5sums
        """
        md5sums = dict()
        self.ftp.retrbinary(f"RETR MD5SUM", open(f"{self.output_dir}/MD5SUM", 'wb').write)
        with open(f"{self.output_dir}/MD5SUM", 'r') as md5sum:
            for line in md5sum:
                hash, file = line.strip().split('  ')
                md5sums[file] = hash
        md5sum.close()
        return md5sums

    def check_md5sum(self, file):
        """
        Check if md5sum of file is matching
        :param file: str, file name
        :return: boolean
        """
        fd = open(f'{self.output_dir}/{file}', 'rb')
        data = fd.read()
        fd.close()
        md5sum = hashlib.md5(data).hexdigest()
        return self.md5sums[file] == md5sum

    def download_epo_alignments(self):
        """
        Download EPO alignments from FTP directory and ensure download was successful by ensuring matching md5sums
        """
        for file in self.files[:5]:# TODO: Remove the [:5]
            complete_path = f"{self.host}/{self.directory}{file}"
            if file == "MD5SUM":
                continue
            elif file.startswith('README'):
                self.ftp.retrbinary(f"RETR {file}", open(f"{self.output_dir}/{file}", 'wb').write)
            # actual EPO alignment files
            else:
                # File exists and md5sum is okay
                if os.path.isfile(f"{self.output_dir}/{file}") and not self.overwrite and self.check_md5sum(file):
                    logging.info(f'{complete_path} exists and md5sum is okay. Skipping.')
                elif os.path.isfile(f"{self.output_dir}/{file}") and not self.overwrite and not self.check_md5sum(file):
                    logging.info(f'{complete_path} exists but md5sum is not matching. Re-downloading the file.')
                elif os.path.isfile(f"{self.output_dir}/{file}") and self.overwrite and self.check_md5sum(file):
                    logging.info(f'Overwriting f"{self.output_dir}/{file}.')
                else:
                    logging.info(f'Downloading {complete_path} to {self.output_dir}/{file}.')
                self.ftp.retrbinary(f"RETR {file}", open(f"{self.output_dir}/{file}", 'wb').write)
                while not self.check_md5sum(file):
                    self.ftp.retrbinary(f"RETR {file}", open(f"{self.output_dir}/{file}", 'wb').write)
                logging.info(f"Successfully downloaded {self.output_dir}/{file}.")


def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('--ftp-dir',
                        help='FTP directory from where to retrieve EPO alignments. '
                             'Default is the 10_primates EPO alignment directory of the current ENSEMBL release '
                             '[http://ftp.ensembl.org/pub/current_emf/ensembl-compara/multiple_alignments/10_primates.epo/]',
                        default='http://ftp.ensembl.org/pub/current_emf/ensembl-compara/multiple_alignments/10_primates.epo/')
    parser.add_argument('--local-epo-alignments-dir',
                        help='Directory to download EPO alignments. Default is the last element of the --ftp-dir, i.e.,'
                             ' 10_primates_epo for the default_parameters', default=None, required=False)
    parser.add_argument('--hg-build',
                        help='The human genome build on which the EPO alignments are based. Default=[GRCh38]',
                        default='GRCh38')
    parser.add_argument('--force', help='Force overwrite of files', action='store_true', default=False)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    epo = EPO(args.ftp_dir, args.local_epo_alignments_dir, args.force)
    epo.download_epo_alignments()


if __name__ == '__main__':
    main(sys.argv[1:])
