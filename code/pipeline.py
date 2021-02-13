#!/usr/lib/python3.8
import argparse, os
from util import cmd

def main():

    parser = argparse.ArgumentParser(description='A pipeline to assembly genomes, recover and annotate variants.')

    subparse = parser.add_subparsers(title='Running options',
                                       description='Choose between running full pipeline or just one specific step.',
                                    dest = 'command' )

    ##### Continuous analysis

    continuous = subparse.add_parser('continuos',
                                       help = 'Will run all sequencial steps in the analysis unless an error or other intereference occurs.',
                                       usage="yet to write")

    continuous.add_argument('--fp','--file_paramet',
                             dest='parameter',
                             type=str,
                             action="store",
                             required=True,
                             default="util/parameters.txt",
                             help="This file is available at <code/util/parameters.txt>, creat your own copy to edit any time by using < cp code/util/parameters.txt ./ > . Use this file to inform the path for files and changes in paramethers. For BWA alignment is mandatory to inform the < path/to/fasta.fa > and < path/to/fastq_folder >. FreeBayes has no mandatory parameters, but you can set target parameter. No mandatory parameter for SnpEff either, default genome reference is the human genome (GRCh37.75).")
     
    

    ##### Discontinuous analysis
    
    discontinuous = subparse.add_parser('discontinuoues',
                                       help = 'Will run separatley the step of you choice BWA, FreeBayes and SnpEff.',
                                       usage='yet to write')

    discontinuous.add_argument('-B','--bwa',
                               dest='bwa',
                               action="store_true",
                               required=False,
                               help="choose running BWA aligner indexer,aligner or both")

    discontinuous.add_argument('-F','--FB',
                               dest='FB',
                               action="store_true",
                               required=False,
                               help="choose running FreeBayes")

    
    discontinuous.add_argument('-S','--SPF',
                                dest='SPF',
                                action="store_true",
                                required=False,
                                help="choose running SnpEff")

    discontinuous.add_argument('--fp','--file_paramet',
                             dest='parameter',
                             type=str,
                             action="store",
                             required=True,
                               default="util/paramet.txt",
                             help=" use this file to inform the path for files and changes in paramethers. For BWA alignment is mandatory to inform path/to/fasta.fa and path/to/fast_folder.") 
     
    discontinuous.add_argument('-s','--sam_file',
                                dest='bwa_output',
                                type=str,
                                action="store",
                                required=False,
                                help="Informa path to <file.sam>. Mandatory for FreeBayes.")

    discontinuous.add_argument('-v','--vcf_file',
                                dest='variant_out',
                                type=str,
                                action="store",
                                required=False,
                                help="Informa path to <file.vcf>. Mandatory for SnpEff, default reference is human genome (GRCh37.75).")

     
    ### RUN PIPELINE
    ###############################################################

    args = parser.parse_args()

    paramt_file = args.parameter

    if hasattr(args,'bwa') or hasattr(args,'FB') or hasattr(args,'SPF'):
        print ('discontinuous')

        bwa = args.bwa
        fb  = args.FB
        spf = args.SPF

        if bwa is True:
            cmd.run_bwa(paramt_file)
            
        if fb is True:
            sam = args.bwa_output
            cmd.run_freebayes(sam,paramt_file)
            
        if spf is True:
            vcf = args.variant_out
            cmd.run_snpeff(vcf,paramt_file)

    else:
        print ('continuous')

        sam = cmd.run_bwa(paramt_file)
        print ('\n\nBWA aligner finished!\n\n')
        vcf = cmd.run_freebayes(sam,paramt_file)
        print ('\n\nFreeBayes variants call finished!\n\n')
        cmd.run_snpeff(vcf, paramt_file)
        print ('\n\nSnpEff functional annotation is finished\n\n')
        print ('Pipeline was a success!"')

        
if __name__ == "__main__":
    main()
                     
