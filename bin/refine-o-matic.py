#!/usr/bin/env python

import sys
import os
import argparse
import logging
import shutil
import itertools
import subprocess
import textwrap
from Bio import SeqIO


class RefineBins:

	FILT_BAM = '_filtered.bam'
	BAM = '.bam'
	FINISHED = '.finished.fasta'
	REASSEMBLED = '.reassembled.fasta'
	def __init__(self):
		self.log_file = 'log.txt'

	def setup(self, output_directory, bin_paths, force):
		print "######### Setting up"

		self.output_directory = output_directory
		try:
			os.mkdir(output_directory)
		except:
			if force:
				shutil.rmtree(output_directory)	
				os.mkdir(output_directory)
			else:
				raise Exception('File exists: %s' % (output_directory))
		self.bin_directory = os.path.join(output_directory,'bin')
		self.proc_directory = os.path.join(output_directory,'proc')
		os.mkdir(self.bin_directory)
		os.mkdir(self.proc_directory)	

		self.log_path = os.path.join(self.output_directory, self.log_file)
		logging.basicConfig(filename=self.log_path, 
							level=logging.INFO, 
							format='%(asctime)s %(levelname)s: %(message)s', 
							datefmt='%m/%d/%Y %I:%M:%S %p')
		logging.info("DIR: %s" % os.getcwd())
		logging.info("CMD: %s" % ' '.join(sys.argv))

		self.contig_to_bin = {}
		for bin in bin_paths:
			cmd = "grep '>' %s | sed 's/>//g' | cut -d' ' -f1 | sort | uniq" \
						% (bin)
			contig_list = set(subprocess.check_output(cmd, shell=True).strip()\
																      .split('\n'))
			for contig in contig_list:
				self.contig_to_bin[contig] = bin

		self.concat_bins = os.path.join(self.output_directory, 'concatenated_bins.fa')
		cmd = 'cat %s > %s' % (' '.join(bin_paths), self.concat_bins)
		logging.info(cmd)
		subprocess.call(cmd, shell=True)

	def _parse_sam_iter(self, sam):
		output_dict = {}
		for line in sam:
			sline = line.strip().split()
			contig_name = sline[2]
			read_name = sline[0]
			bin_id = self.contig_to_bin[contig_name]
			yield bin_id, read_name

	def get_bam(self, directory):
		bam_list 	  = [] 
		filt_bam_list = []
		for file in os.listdir(directory):
			if file.endswith(self.BAM):
				if file.endswith(self.FILT_BAM):
					filt_bam_list.append(file)
				else:
					bam_list.append(file)

		if any(filt_bam_list):
			return filt_bam_list
		else:
			return bam_list

	def collect_bin_stats(self, bin):
		contigs = list(SeqIO.parse(open(bin), 'fasta'))
		ambiguous_bases = 0
		total_bases     = 0
		average_list 	= [] 
		n_seq 			= len(contigs)
		for record in contigs:
			sequence = str(record.seq)
			total_bases+=len(sequence)
			average_list.append(len(sequence))
			ambiguous_bases+=sequence.count('N')

		if bin.endswith(self.REASSEMBLED):
			bin_type = 'reassembled'
		elif bin.endswith(self.FINISHED):
			bin_type = 'finished'
		else:
			bin_type = 'input'

		return [bin_type,
				str(n_seq), 
				str(total_bases), 
				str(ambiguous_bases),
				str(min(average_list)), 
				str(max(average_list)), 
				str(total_bases/float(n_seq))]

	def mapping(self, forward_reads, reverse_reads, threads, bam_filter=None):
		print "######### Mapping"
		cmd = 'bamm make --quiet -o %s -k -d %s -c %s -t %s 2>/dev/null' \
			% (self.output_directory, 
			   self.concat_bins, 
			   ' '.join(list(itertools.chain(*zip(forward_reads,
			   							 		  reverse_reads)))),
			   threads)

		logging.info(cmd)
		subprocess.call(cmd, shell=True)

		if bam_filter:

			cmd = "bamm filter -b %s -o %s --percent_id %s" % (self.get_bam(self.output_directory),
															   self.output_directory,
														 	   bam_filter)
			logging.info(cmd)
			subprocess.call(cmd, shell=True)

		bin_to_reads = {}
		for bam in self.get_bam(self.output_directory):
			cmd = 'samtools view %s' % os.path.join(self.output_directory, bam)
			for bin, read in self._parse_sam_iter(subprocess.check_output(cmd, shell=True).strip().split('\n')):
				if bin not in bin_to_reads:
					bin_to_reads[bin] = [read]
				else:
					bin_to_reads[bin].append(read)

		mapped_read_names = list(itertools.chain(*bin_to_reads.values()))

		# Extract mapped reads out of original reads file
		output_forward_reads = os.path.join(self.output_directory, 'all_mapped.1.fastq')
		output_reverse_reads = os.path.join(self.output_directory, 'all_mapped.2.fastq')
		
		for for_pair, rev_pair in zip(forward_reads, reverse_reads):
			cmd = 'fxtract -HXf /dev/stdin %s >> %s ' % (for_pair, output_forward_reads) 
			logging.info(cmd)
			process = subprocess.Popen(["bash", "-c", cmd], 
										stdin=subprocess.PIPE,
										stdout=subprocess.PIPE)
			process.communicate('\n'.join(mapped_read_names))
			cmd = 'fxtract -HXf /dev/stdin %s >> %s ' % (rev_pair, output_reverse_reads)
			logging.info(cmd)
			process = subprocess.Popen(["bash", "-c", cmd], 
										stdin=subprocess.PIPE,
										stdout=subprocess.PIPE)
			process.communicate('\n'.join(mapped_read_names))
		
		for bin_path, reads_list in bin_to_reads.items():
			bin_basename = os.path.splitext(os.path.basename(bin_path))[0]
			bin_forward_reads = os.path.join(self.output_directory,
											 bin_basename + '_mapped.1.fastq')
			bin_reverse_reads = os.path.join(self.output_directory,
											 bin_basename + '_mapped.2.fastq')
			cmd = 'fxtract -HXf /dev/stdin %s > %s' \
							% (output_forward_reads, bin_forward_reads)
			logging.info(cmd)

			process = subprocess.Popen(["bash", "-c", cmd], 
										stdin=subprocess.PIPE,
										stdout=subprocess.PIPE)
			process.communicate('\n'.join(reads_list))
			cmd = 'fxtract -HXf /dev/stdin %s > %s' \
							% (output_reverse_reads, bin_reverse_reads)
			logging.info(cmd)
			process = subprocess.Popen(["bash", "-c", cmd], 
										stdin=subprocess.PIPE,
										stdout=subprocess.PIPE)
			process.communicate('\n'.join(reads_list))


			yield bin_forward_reads, bin_reverse_reads, bin_path, bin_basename

	def reassembly(self, forward, reverse, bin_reference_contigs,
					bin_name, threads):
		print "######### Re-assembly"

		assembly_output_directory = os.path.join(self.output_directory, bin_name + '.spades')
		cmd = 'spades.py -o %s -t %s -1 %s -2 %s --trusted-contigs %s --careful --phred-offset 33 1> /dev/null' \
						% (assembly_output_directory, threads, forward, reverse, bin_reference_contigs)

		logging.info(cmd)
		subprocess.call(cmd, shell=True)

		reassembled_bin_path = os.path.join(self.bin_directory, '%s.reassembled.fasta' % bin_name)
		shutil.copy(os.path.join(assembly_output_directory, 'contigs.fasta'), 
					reassembled_bin_path)
		return reassembled_bin_path

	def finishm(self, bin_contigs, bin_name, forward_reads, reverse_reads):
		print "######### Running FinishM" 
		finishm_output_directory = os.path.join(self.output_directory, 'finishM')
		velvet_output_directory = os.path.join(finishm_output_directory, "velvet")
		cmd = 'finishm roundup --quiet --genomes %s --fastq %s,%s --output-directory %s --velvet-directory %s' \
					% (bin_contigs, forward_reads, reverse_reads, finishm_output_directory, velvet_output_directory)
		logging.info(cmd)
		subprocess.call(cmd, shell=True)
		finished_bin = os.path.join(finishm_output_directory, '%s.reassembled.fasta.scaffolds.fasta' % bin_name)
		finished_renamed_bin = os.path.join(self.bin_directory, '%s.finished.fasta' % bin_name)
		shutil.copy(finished_bin, finished_renamed_bin)
		return finished_renamed_bin

	def write_stats(self, output_file, stats_lines):
		header = ['bin', 'bin_type', 'number_of_sequences', 'total_bases', 'ambiguous_bases',
				  'smallest_contig', 'largest_contig', 'average_contig_size']
		with open(output_file, 'w') as out_io:
			out_io.write('\t'.join(header) +'\n')
			for line in stats_lines:
				out_io.write('\t'.join(line) +'\n')

	def cleanup(self):
		for file in os.listdir(self.output_directory):
			if(file.startswith('concatenated_bins.') or
			   file.startswith('all_mapped') or
			   file.endswith(self.BAM) or
			   file.endswith('.fastq')):
				source = os.path.join(self.output_directory, file)
				destination = os.path.join(self.proc_directory, file)
				shutil.move(source, destination)
	def run_checkm(self, threads):
		print "######### Running checkM on bins"

		checkm_output_directory = os.path.join(self.output_directory, 'checkm_lineage_wf')
		checkm_output_file = os.path.join(self.output_directory, 'checkm_results.tsv')
		
		cmd = 'checkm lineage_wf -q -t %s -x fasta --tab_table -f %s %s %s' \
				% (threads, checkm_output_file, self.bin_directory, checkm_output_directory)

		logging.info(cmd)
		subprocess.call(cmd, shell=True)

	def main(self, args):
		print "\n######### Starting refine-o-matic"

		stats_lines = []
		# Setup output directory
		self.setup(args.output_directory, args.bins, args.force)

		# mapping
		for extracted_forward_reads, extracted_reverse_reads, bin_path, bin_name in self.mapping(args.forward, 
																								 args.reverse,
																								 args.threads,
																								 args.bam_filter):
			print "######### Working on %s" % bin_name		
			# spades.py re-assembly
			reassembled_bin_path = self.reassembly(extracted_forward_reads,
												   extracted_reverse_reads, 
									 			   bin_path,
									  			   bin_name,
								 	  			   args.threads)

			# finishM
			finished_bin_path = self.finishm(reassembled_bin_path,
											 bin_name, 
											 extracted_forward_reads, 
											 extracted_reverse_reads)

			stats_lines.append([os.path.basename(bin_path)] + self.collect_bin_stats(bin_path))
			stats_lines.append([os.path.basename(reassembled_bin_path)] + self.collect_bin_stats(reassembled_bin_path))
			stats_lines.append([os.path.basename(finished_bin_path)] + self.collect_bin_stats(finished_bin_path))

			print "######### Finished %s" % bin_name
		
		print "######### Writing simple bin stats" 
		self.write_stats(os.path.join(self.output_directory, 'stats.tsv'), stats_lines)
		
		if args.checkm:
			self.run_checkm(args.threads)
		
		self.cleanup()
		print "######### refine-o-matic has finished. Beep-bop-boop.\n"

if __name__=="__main__":
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
									 description=textwrap.dedent('''
			Refine-o-matic: Just a simple pipeline for 
			lazy people who want refine their bins

			Pipeline
			----------------------------------------
			1) map reads to bins (BamM)
			2) re-assemble mapped reads (Spades)
			3) Gap filling and scaffolding (FinishM)'''))

	parser.add_argument('--bins', nargs='+',  help='space separated list of bins to refine', required=True)
	parser.add_argument('--forward', nargs='+', help='forward reads', required=True)
	parser.add_argument('--reverse', nargs='+', help='reverse reads', required=True)
	parser.add_argument('--output_directory', help='output directory', required=True)
	parser.add_argument('--threads', default='10', help='number of threads to use')
	parser.add_argument('--force', help='overwrite previous run', action='store_true')
	parser.add_argument('--checkm', help='run checkM on refined bins', action='store_true')
	parser.add_argument('--bam_filter', help='filter bam by percent id mapping', type = float)

	args = parser.parse_args()

	rb = RefineBins()
	rb.main(args)
