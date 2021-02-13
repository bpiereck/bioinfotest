import os


################# accessing parameters

def bwa_paramet(file):

    PARAMETERS = {}                 # Parameters dictionary
    
    PARAMETERS['thread'] = ' '      # default No threads config
    PARAMETERS['tipo']   = 'bwtsw'  # default big genomes
    PARAMETERS['create'] = 'Yes'    # default creat index before running
            
    with open(file,'r') as f:
        for line in f:
            line = line.strip()

            if line.startswith('bwa'):                            # check for bwa parameters only
                header, paramt = line.split('=')
                paramt = paramt.strip()
                if 'fasta' in line:                               # get path to fast
                    PARAMETERS['fasta'] = paramt.strip()
                    if '.fasta' in paramt:
                      paramt = paramt.strip()
                      prefix_df = paramt[:-6]                     # default prefix if file.fasta
                      PARAMETERS['prefix'] = prefix_df
                    elif '.fa' in paramt:
                      paramt = paramt.strip()
                      prefix_df = paramt[:-3]                     # default prefix if file.fa
                      PARAMETERS['prefix'] = prefix_df
                elif 'idxType' in line:
                    PARAMETERS['tipo'] = paramt.strip()           # get index type bwtsw/is
                elif 'output' in line:
                    PARAMETERS['output'] = paramt.strip()         # bwa output name
                elif 'fastq' in line:
                    PARAMETERS['fastq'] = paramt.strip()          # get path for fastq-containing folder
                elif 'thread' in line:
                    PARAMETERS['thread'] = f'-t {paramt.strip()}' # save nuber of threads if param is given
                elif 'prefix' in line:
                    PARAMETERS['prefix'] = paramt.strip()         # save prefix name if given
                elif 'create_' in line:
                    PARAMETERS['create'] = paramt.strip()         # create idex option Yes/yes/No/no
                

    return PARAMETERS


def fb_paramet(file):

    PARAMETERS = {}                          # Parameters dictionary
    PARAMETERS['target'] = ' '               # default no target, therefore whole genome
        
    with open(file,'r') as f:
        for line in f:
            line = line.strip()

            if line.startswith('fb'):                           # check for FreeBayes parameters only
                header, paramt = line.split('=')
                if 'target' in line:
                    PARAMETERS['target'] = f'--target {paramt.strip()}' # save path to target if given
                if 'fasta' in line:
                    PARAMETERS['fasta'] = paramt.strip()                # get path to fast

    return PARAMETERS


def snp_paramet(file):

    PARAMETERS = {}                      # Parameters dictionary
    PARAMETERS['thread'] = ' '           # default No threads config
    PARAMETERS['genome'] = 'GRCh37.75'   # default human genome
    PARAMETERS['mem'] = '-Xmx8g'         # default method
    
    with open(file,'r') as f:
        for line in f:
            line = line.strip()

            if line.startswith('snp'):                     # check for SnpEff parameters only
                header, paramt = line.split('=')
                if 'memmory' in line:                      
                    PARAMETERS['mem'] = paramt.strip()             # get memmory use limit
                elif 'genome' in line:
                    PARAMETERS['genome'] = paramt.strip()          # get ref. genome if given
                elif 'thread' in line:
                    PARAMETERS['thread'] = f'-t {paramt.strip()}'  # get thread if given
                    

    return PARAMETERS

###################################################3


def run_bwa(paramt):
    
    PARAMETERS  = bwa_paramet(paramt)    # Get parameters dictionary
    prefix      = PARAMETERS['prefix']   # prefix
    tipo        = PARAMETERS['tipo']     # running type for indexing
    fasta       = PARAMETERS['fasta']    # fasta file
    bwa_output  = PARAMETERS['output']   # output name file.sam
    path_fastq  = PARAMETERS['fastq']    # path to fastq folder
    thread      = PARAMETERS['thread']   # running threads
    create      = PARAMETERS['create']   # Yes/yes if index does't exist OR No/no if index already exist

    print (f'''\n\n
##################################
\BWA LIST OF PARAMETERS:
##################################
threads      = {thread}
create index = {create}
Index prefix = {prefix}
bwa run type = {tipo}
fastq folder = {path_fastq}

-- input  file  = {fasta}
-- output file  = {bwa_output}\n\n
''')

    
    if create == 'Yes' or create == 'yes':                     # if index does't exist creat one
      
      cmd_index  = f'bwa index -p {prefix} -a {tipo} {fasta}'  # creat index cmd
      print (f'bwa build index cmd:\n{cmd_index}\n\n')
      os.system(cmd_index)                                     # run bwa index option
      

    if create == 'No' or create == 'no' or create == 'Yes' or create == 'yes':  # if index exists, save time running directly the assembly
        
        listdir = os.listdir(path_fastq)                           # list all files in given path
        fastq   = []                                               # make list of fastq files
        for file in listdir:                                       
          if file.endswith('fastq') or file.endswith('fastq.gz'):  # get only fastq files
              file = os.path.join(path_fastq,file)                 # join path/file.fastq
              fastq.append(file)                                   # add to list
              
        fastq = ' '.join(fastq)                                         # use list to creat string of fastq files
        cmd_alig  = f'bwa mem {thread} {prefix} {fastq} > {bwa_output}' # cmd to run BWA assembly
        print (f'bwa assembly cmd:\n{cmd_alig}\n\n')
        os.system(cmd_alig)                                             # run BWA assembly
    

    return bwa_output
    



def run_freebayes(bwa_output,paramt):

    PARAMETERS  = fb_paramet(paramt)      # Get parameters dictionary
    target      = PARAMETERS['target']    # target file
    fasta       = PARAMETERS['fasta']     # fasta file
    nameB       = bwa_output[:-4]
    bam_file    = f'{nameB}.bam'          # create bam file name
    nameV       = bam_file[:-4]
    variant_out = f'{nameV}_var.vcf'      # creat vcf output name

    print (f'''\n\n
##################################
FreeBayes LIST OF PARAMETERS:
##################################
target       = {target}
input  file  = {fasta}
output bam   = {bam_file}
output vcf   = {variant_out}\n\n
''')
    
    if os.path.exists(bwa_output) is True:                                # if bwa.sam file was correctly generated
        
        cmd_sam_bam = (f'samtools view -S -b {bwa_output} > {bam_file}')  # cmd to convert SAM to BAM file
        print (f'convert SAM to BAM cmd:\n{cmd_sam_bam}\n\n')
        os.system(cmd_sam_bam)                                            # run Samtool to make conversion

        # sort and index bam file  to use in Freebaye
        cmd_bam_index = f'samtools sort {bam_file} -o {bam_file}.sort ; mv {bam_file}.sort {bam_file}; samtools index {bam_file}' 
        
        print (f'Sort and index BAM cmd:{cmd_bam_index}\n\n')
        os.system(cmd_bam_index)       # run Samtool to creant BAM index -> file.bai
        bai_file = f'{bam_file}.bai'   # name bai file to run next step

        
        if os.path.exists(bam_file) is True and os.path.exists(bai_file) is True:             # if both files are correctly generated

            cmd_freebayes = (f'freebayes -f {fasta} {bam_file} {target} > {variant_out}')  # cmd to call variantes (vcf file)
            print (f'Call variantes cmg:\n{cmd_freebayes}')
            os.system(cmd_freebayes)                                                          # run freebayes and get variantes
            
        else:
            if os.path.exists(bam_file) is False:
                    print ('ERROR could not find {bam_file}')
            elif os.path.exists(bai_file) is False:
                    print ('ERROR could not find {bai_file}')
        
    else:
        print (f'ERROR could not find assembly file: {bwa_output}')

    return variant_out




def run_snpeff(variant_out,paramt):

    PARAMETERS = snp_paramet(paramt)
    thread    = PARAMETERS['thread']  # get number of threads if given, otherwise none
    gen_ref    = PARAMETERS['genome'] # get reference, otherwise homo sapiens
    mem        = PARAMETERS['mem']    # default mem limit -Xmx8g

    print (f'''\n\n
##################################
snpEff LIST OF PARAMETERS:
##################################
threads   = {thread}
reference = {gen_ref} genome
memmory   = {mem}
output    = {variant_out}.annotation
''')
    
    ### adicionar
    # snpEff download -v GRCh37.75 # download databse if not available
    # snpEff -Xmx8g -c ./snpEff.config -v GRCh37.75 # cinfig databse before using


    if os.path.exists(variant_out) is True:
          cmd_snpeff = f'snpEff {mem} -v {gen_ref} {variant_out} > {variant_out}.annotation'
          print (f'Annotation of variantes cmd:\n{cmd_snpeff}\n\n')
          os.system(cmd_snpeff)
    else:
        (f'ERROR could not find vcf file: {variant_out}')
    



