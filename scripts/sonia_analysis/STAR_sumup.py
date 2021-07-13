import os,sys
import pandas as pd

def read_log_final_out(file,name):
    log_final_out_df=pd.read_csv(file,sep='|',header=None, names=['info',name])
    log_final_out_df['info']=log_final_out_df['info'].str.strip()
    log_final_out_df[name]=log_final_out_df[name].str.strip()
    return log_final_out_df


def main():
    folder=sys.argv[1]
    output=sys.argv[2]
    os.chdir(folder)
    command='find -name "*.final.out" -print0 | xargs -0 ls > file_list'
    os.system(command)
    with open('file_list','r') as f:
        for i,line in enumerate(f):
            file=line.rstrip('\n')
            name=file.split('/')[-2]
            log_final_out_df=read_log_final_out(file,name)
            if i==0:
                log_final_out_df_final=log_final_out_df.copy(deep=True)
            else:
                log_final_out_df_final=pd.merge(log_final_out_df_final,log_final_out_df,how='left',on='info')
    os.system('rm file_list')
    log_final_out_df_final.to_csv(path_or_buf=output,
                  index=False, header=True, sep=',')



if __name__=='__main__':
    main()
