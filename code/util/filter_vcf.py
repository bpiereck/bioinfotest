

def get_impact(info):

    MUTATION = []
    
    info   = info.split('|')
    size   = len(info)
    ciclos = size/15
    cicle_pos = 2
    i = 0
    while i < ciclos-1:
        n = int(cicle_pos-1)
        mutation = info[n].strip()
        impact   = info[cicle_pos].strip()
        if impact != "MODIFIER":
          MUTATION.append((mutation, impact))
        cicle_pos += 15
        i += 1

    return MUTATION


import sys


def filtering(vcf):

  filtered = {} # index: subdict
  subdict  = {} # 'chr':chr, 'posix':posix, 'mutations':[(mut,imp),(mut,imp)]
  index = 1
  
  with open(vcf,'r') as vcf:
      
    for line in vcf:
      line = line.strip()
      
      if 'ANN=' in line:
        if 'HIGH' in line or 'LOW' in line or 'MODERATE' in line:

          subdict  = {} #'chr':chr,'posix':posix,'mutations':[(m,i),(m,i)]  
          columns = line.split('\t')
          
          Chr = columns[0]
          subdict['chr'] = Chr
          
          posix = columns[1]
          subdict['posix'] = posix
          
          info  = columns[7]
          mutations = get_impact(info)
          subdict['mutations'] = mutations
          
          filtered[index] = subdict
          index += 1
          

    return filtered

def filter(vcf):
    
    filtered=filtering(vcf)

    with open('ANSWERS.annot.txt','w') as out:

      str = f'SEQ_ID\tCHR\tPOSIÇÃO\tMUTAÇÕES\tIMPACTO\n'
      out.write(str)
        
      for index, value in filtered.items():
        Chr   = value['chr']
        posix = value['posix']
        mutations = value['mutations']


        for tup in mutations:
          mutation = tup[0]
          impact   = tup[1]

          location = f'{index}\t{Chr}\t{posix}\t{mutation}\t{impact}\n'
          out.write(location)
          


        

    
vcf = sys.argv[1]
filter(vcf)
