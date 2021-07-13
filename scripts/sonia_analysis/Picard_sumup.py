import os,sys
import pandas as pd

pd.set_option('max_rows', 100)
pd.set_option('max_columns', 80)
pd.set_option('display.width', 400)

def read_Picard_file(file,name):
    with open(file,'r') as p:
        Picard_content=p.readlines()
    start_index=[i for i in range(len(Picard_content)) if Picard_content[i].startswith('MEDIAN_INSERT_SIZE')==True][0]
    columns=list(Picard_content[start_index].rstrip().split())[:6]
    values=list(Picard_content[start_index+1].rstrip().split())[:6]
    Picard_df=pd.DataFrame(data={'info':columns,name:values})
    return Picard_df


def main(folder='',output=''):
    if folder=='':
        folder=sys.argv[1]
        output=sys.argv[2]
    os.chdir(folder)
    command='find -name "*Picard_insert_size_metrics.txt" -print0 | xargs -0 ls > file_list'
    os.system(command)
    with open('file_list','r') as f:
        for i,line in enumerate(f):
            file=line.rstrip('\n')
            name=file.split('/')[-2]
            Picard_df=read_Picard_file(file,name)
            if i==0:
                Picard_df_final=Picard_df.copy(deep=True)
            else:
                Picard_df_final=pd.merge(Picard_df_final,Picard_df,how='left',on='info')

    os.system('rm file_list')

    cols = list(Picard_df_final)
    cols.insert(0, cols.pop(cols.index('info')))
    Picard_df_final = Picard_df_final.loc[:, cols]
    print(Picard_df_final[:5])

    Picard_df_final.to_csv(path_or_buf=output,
                  index=False, header=True, sep=',')



if __name__=='__main__':
    folder=sys.argv[1]
    output=sys.argv[2]
    #folder='/home/vincent/Desktop/Sequencing/Tissues/Picard'
    #output='/home/vincent/Desktop/Sequencing/Tissues/Picard/Picard_sumup.csv'
    main(folder,output)
