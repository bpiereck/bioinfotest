import sys, re, os
import pickle #implement binary protocol

### Goal: montar um cromossomo teorico a partir do alinhamento de sequências
## string de até 1000 bases
## sobreposição mais da metade
## menor contig resultante com todas as strings



def make_index(db):

    print(f'\n\nmaking index of <{db}> .............')

    Index = {} # cread dict as --> seq:{index:(start,end,idx)}
    
    with open(db, 'r') as f:
        
        for index, seq in enumerate(f,1):    # enumerate lines start in 1
            seq = seq.strip()                # clean edges

            if seq not in Index:
                if not seq:
                    pass
                else:
                    #start,end = frag_third(seq)   # get 2/3 and last 1/3
                    #Index[seq]= (start,end,index)  # add to dict
                    Index[seq]= index
                                
        pickle.dump(Index, open(db + '.index', "wb")) # save it binary
        print (f'DONE!\n')
        
    

def alignment(db):

    index= f'{db}.index'                   # make index name
    INDEX = pickle.load(open(index,'rb'))  # load and read index binary

    print(f'starting alignmet of <{db}>')
    
    SEQ = ""                               # build big SEQ
    assembled = []                         # save idx of assembled seq so is not read twice
    

    while len(assembled) < len(INDEX):     # try to find match untill all reads are part of contig
      for key, value in INDEX.items():
        seq   = key                        # get seq to align
        idx   = value                      # idx to control assembled reads

        if idx not in assembled:           # if not yet assembled
          if SEQ == '':                    # if it is the 1st read
            SEQ = seq                      # make it main SEQ
            assembled.append(idx)          # save as already assembled
            
                                  
          else:                            # for all other try the following
            third = int(len(seq)/3)
            N = third*2
            half = int(len(seq)/2)
                    
            third = int(len(seq)/3)
            N = third*2
            if seq[:N] in SEQ:             # if initial 2/3s matchs with SEQ
              start = seq[:N]              # get initial 2/3rds of seq to alig 
              end   = seq[N:]              # save final 1/3rd
              result = re.search(rf'{start}',SEQ, re.IGNORECASE)  # find alignment 
              a,b = (result.span())                               # get alignment location
              SEQ = f'{SEQ[0:a]}{start}{end}'                     # concatenate both reads              
              assembled.append(idx)                               # save as assembled
              

            elif seq[third:] in SEQ:   # if final 2/3s matchs with SEQ
              start = seq[:third]      # get initial 1/3rd of seq
              end   = seq[third:]      # get final 2/3rds of seq to alig
              result = re.search(rf'{end}',SEQ, re.IGNORECASE)  # find alignment 
              a,b = (result.span())                             # get alignment location
              SEQ = f'{start}{end}{SEQ[b:]}'                    # concatenate both reads              
              assembled.append(idx)                             # save as assembled
              

            else:                      # if can't match 2/3rds, try 1/2 seq
              half = int(len(seq)/2)
              if seq[:half+1] in SEQ:  # if inital 1/2 in matchs with SEQ
                start = seq[:half+1]   # get initial half, assure if odd number 1st half is bigger
                end   = seq[half+1:]   # get final half
                result = re.search(rf'{start}',SEQ, re.IGNORECASE)  # find alignment
                a,b = (result.span())                               # get alignmet location
                SEQ = f'{SEQ[0:a]}{start}{end}'                   # concatenate both reads
                assembled.append(idx)                               # save as assembled
                
                  
              elif seq[half-1:] in SEQ:      # if final1/2 matches with SEQ
                start = seq[:half-1]         # get initial half
                end   = seq[half-1:]         # get final half, assure if odd number last half is bigger
                result = re.search(rf'{end}',SEQ, re.IGNORECASE)    # find alignment
                a,b = (result.span())                               # get alignmet location                  
                SEQ = f'{start}{end}{SEQ[b:]}'                      # concatenate both reads
                assembled.append(idx)                               # save as assembled
                

              else:
                half = int(len(seq)/2)
                if seq[:half] in SEQ:    # if inital 1/2 in matchs with SEQ
                    start = seq[:half]   # get initial half, assure if odd number 1st half is bigger
                    end   = seq[half:]   # get final half
                    result = re.search(rf'{start}',SEQ, re.IGNORECASE)  # find alignment
                    a,b = (result.span())                               # get alignmet location
                    SEQ = f'{SEQ[0:a]}{start}{end}'                   # concatenate both reads
                    assembled.append(idx)                               # save as assembled
                    
                  
                elif seq[half:] in SEQ:      # if final1/2 matches with SEQ
                  start = seq[:half]         # get initial half
                  end   = seq[half:]         # get final half, assure if odd number last half is bigger
                  result = re.search(rf'{end}',SEQ, re.IGNORECASE)    # find alignment
                  a,b = (result.span())                               # get alignmet location                  
                  SEQ = f'{start}{end}{SEQ[b:]}'                      # concatenate both reads
                  assembled.append(idx)                               # save as assembled

                  
    out = db.split('.')
    out = '.'.join(out[:-1])
    out = f'{out}.out'
    with open(out,'w') as out:
        out.write(SEQ)

    os.system(f'rm {db}.index')

    print ('alignment DONE!\n\n')
    return SEQ

                
    
            
            
            
                                

                                     
    
    

###############################
###############################
###############################

db = sys.argv[1]
make_index(db)
alignment(db)
