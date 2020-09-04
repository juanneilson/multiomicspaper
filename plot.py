import numpy as np
import matplotlib.pyplot as plt

n_modifications = 3


def plot_distribution_of_designs(df): 
    bar_height = 1
    labels = ['KO', 'NoMod', 'UP']
    colors = ['#019600', 'grey', '#219AD8']
        
    plt.style.use('seaborn-white')
    
    dataframe = df.copy()
    reactions = dataframe.columns
    
    n_rec = len(dataframe)
    dataframe.loc[n_rec] = [[list(dataframe[reaction]).count(int(i))/n_rec*100 
                             for i in range(n_modifications)]  for reaction in reactions]
    
    data = [ [dataframe.iloc[-1][r][num] for r in reactions] 
            for num in range(n_modifications)]
    
    y_pos = np.arange(len(reactions))

    fig = plt.figure(figsize=(7,5))
    ax = fig.add_subplot(111)

    # Remove frame
    for spine in plt.gca().spines.values():
        spine.set_visible(False)

    patch_handles = []
    # left alignment of data starts at zero
    left = np.zeros(len(reactions)) 
    for i, d in enumerate(data):
        patch_handles.append(ax.barh(y_pos, d, 
                                     color=colors[i%len(colors)], edgecolor='white',
                                     height=bar_height, align='center', 
                                     left=left, label=labels[i]))
        left += d

    # search all of the bar segments and annotate
    for j in range(n_modifications):
        for i, patch in enumerate(patch_handles[j].get_children()):
            bl = patch.get_xy()
            x = 0.5*patch.get_width() + bl[0]
            y = 0.5*patch.get_height() + bl[1]
            ax.text(x,y, "%d%%" % (data[j][i]), ha='center')

    ax.set_title('Distribution of modifications')
    plt.tick_params(top='off', bottom='off', left='off', right='off', labelleft='on', 
                    labelbottom='off')
    plt.yticks(y_pos, reactions)
    ax.invert_yaxis()
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.show()
    
    
def plot_DO_extmets(od,ext_metabolites):
    fig, ax = plt.subplots(figsize=(12,4), ncols=2, nrows=1)
    od.plot(ax=ax[0], style='s-', title='Cell', label='dcw', legend=True)
    ax[0].set_xlabel("Hour")
    ax[0].set_ylabel("Concentration [gDW/L]")
    ext_metabolites.plot(ax=ax[1], style='o-', title='External Metabolites')
    ax[1].set_xlabel("Hour")
    ax[1].set_ylabel("Concentration [mM]")
