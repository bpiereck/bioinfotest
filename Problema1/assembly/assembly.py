import sys, re, os
import pickle #implement binary protocol

### Goal: montar um cromossomo teorico a partir do alinhamento de sequências
## string de até 1000 bases
## sobreposição mais da metade
## menor contig resultante com todas as strings


def frag_third(seq):
    third = int(len(seq)/3)
    N = third*2
    start  = seq[0:N]    # get 2/3
    end    = seq[N:]     # get last 1/3

    return start, end


def frag_half(seq):
    half = int(len(seq)/2)
    start  = seq[0:half]  # get 1/2
    end    = seq[half:]   # get last 1/2

    return start, end



def make_index(db):

    Index = {} # cread dict as --> seq:{index:(start,end,idx)}
    
    with open(db, 'r') as f:
        
        for index, seq in enumerate(f,1):    # enumerate lines start in 1
            seq = seq.strip()                # clean edges

            if seq not in Index:
                if not seq:
                    pass
                else:
                    start,end = frag_third(seq)   # get 2/3 and last 1/3
                    Index[seq]=(start,end,index)  # add to dict
                                
        pickle.dump(Index, open(db + '.index', "wb")) # save it binary
        
    

def alignment(db):

    index= f'{db}.index'                   # make index name
    INDEX = pickle.load(open(index,'rb'))  # load and read index binary

    SEQ = ""                               # build big SEQ
    assembled = []                         # save idx of assembled seq so is not read twice
    

    while len(assembled) < len(INDEX):     # try to find match untill all reads are part of contig
      for key, value in INDEX.items():
        seq   = key                        # get seq to align
        start = value[0]                   # first 2/3 of seq to compare
        end   = value[1]                   # last 1/3 of seq
        idx   = int(value[2])              # idx to control assembled reads

        if idx not in assembled:           # if not yet assembled
          if SEQ == '':                    # if it 1st read
            SEQ = seq                      # make it main SEQ
            assembled.append(idx)          # save as already assembled
            
           
                        
          else:
            if start in SEQ:                                          # if 2/3s matchs with SEQ
              result = re.search(rf'{start}',SEQ, re.IGNORECASE)      # find alignment 
              assembled.append(idx)                                   # save as assembled
              a,b = (result.span())                                   # get alignment location
              SEQ = f'{SEQ[0:a+1]}{start}{end}'                       # concatenate both reads
                
             
            elif start not in SEQ:                                    # if can't find 2/3 match, try 1/2 match
                start, end  = frag_half(seq)                          # get new star and end 1/2

                if start in SEQ:                                      # if 1/2 in matchs with SEQ
                  result = re.search(rf'{start}',SEQ, re.IGNORECASE)  # find alignment
                  assembled.append(idx)                               # save as assembled
                  a,b = (result.span())                               # get alignmet location
                  SEQ = f'{SEQ[0:a+1]}{start}{end}'                   # concatenate both reads
                  

    out = db.split('.')
    out = '.'.join(out[:-1])
    out = f'{out}.out'
    with open(out,'w') as out:
        out.write(SEQ)

    os.system(f'rm {db}.index')

    print ('alignment has finished')
    return SEQ

                
    
            
            
            
                                

                                     
    
    

###############################
###############################
###############################

db = sys.argv[1]
make_index(db)
alignment(db)
