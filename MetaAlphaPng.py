import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator, LinearLocator
import pandas as pd

# plt.rcParams['font.sans-serif'] = ['SimHei']


def alph_png(alph,out_name):
    alph=alph.set_index('α多样性')
    cankao=alph['参考值'].astype(float)
    del alph['参考值']
    #alph['T2']=[3.5,0.9]
    #alph['T3']=[3.2,0.6]
    alph=alph.T
    fig, axs = plt.subplots(1, 2, figsize=(4, 1.5),dpi=300)
    num_T=-1
    line_num=-1
    for i in alph['Shannon香农指数']:
        #print(i)
        num_T+=1
        x=num_T
        if i=='':
            axs[0].scatter(x,y,color='white',s=0)
        else :
            line_num+=1
            y=i
            axs[0].scatter(x,y,color='#0AEE1A',s=8)
            if line_num>0:
                axs[0].plot([x1,x],[y1,y],color='#0AEE1A',linewidth=0.5)
            #用于连线
            x1=num_T
            y1=i
    num_T=-1
    line_num=-1
    for ia in alph['Simpson辛普森指数']:
        #print(i)
        num_T+=1
        x=num_T
        if ia=='':
            axs[1].scatter(x,y,color='white',s=0)
        else :
            line_num+=1
            y=ia
            axs[1].scatter(x,y,color='#24A4EA',s=8)
            if line_num>0:
                axs[1].plot([x1,x],[y1,y],color='#24A4EA',linewidth=0.5)
            #用于连线
            x1=num_T
            y1=ia
    #去掉边框
    axs[0].spines['right'].set_visible(False)
    axs[0].spines['top'].set_visible(False)
    axs[1].spines['right'].set_visible(False)
    axs[1].spines['top'].set_visible(False)
    #刻度数量
    axs[0].yaxis.set_major_locator(LinearLocator(3))
    axs[1].yaxis.set_major_locator(LinearLocator(3)) 
    #坐标轴刻度字体大小
    axs[0].tick_params(labelsize=8)
    axs[1].tick_params(labelsize=8)
    axs[0].set_xticks(range(len(alph['Shannon香农指数'].index)))
    axs[0].set_xticklabels(alph['Shannon香农指数'].index.tolist())
    axs[1].set_xticks(range(len(alph['Simpson辛普森指数'].index)))
    axs[1].set_xticklabels(alph['Simpson辛普森指数'].index.tolist())
    # axs[1].set_xticklabels(['','T1','T2','T3'])
    #标题
    axs[0].set_title('Shannon',size=8)
    axs[1].set_title('Simpson',size=8)
    #参考值线
    axs[0].plot((0,2),(cankao.loc['Shannon香农指数'],cankao.loc['Shannon香农指数']),linewidth=0.7,color='black',linestyle='-.')
    axs[1].plot((0,2),(cankao.loc['Simpson辛普森指数'],cankao.loc['Simpson辛普森指数']),linewidth=0.7,color='black',linestyle='-.')
    #子图间隔
    fig.subplots_adjust( wspace=0.4)
    fig.savefig(out_name,transparent=True,bbox_inches='tight')
    # plt.subplots_adjust( wspace=0.4)
    # plt.savefig(out_name,transparent=True,bbox_inches='tight')
    # plt.close('all')