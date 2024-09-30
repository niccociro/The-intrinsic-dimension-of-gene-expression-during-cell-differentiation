import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np



grid_linewidth = 0.5
grid_linestyle = '--'
grid_axis = 'both'


def plot_temporal_trend(results_dict, dataset_name, xlabel,
               x_ordering = None, title = '',
               plot_params = {"figsize": (15, 10),
                              "label_size": 30,
                              "markersize": 20,
                              "linewidth": 2,
                              "xlabel_rotation": 0}):
    
        color = 'maroon'

        x = results_dict['Mean trend']['Labels']
        y = results_dict['Mean trend']['Mean']
        y_err = results_dict['Mean trend']['Std']

        if(x_ordering is not None):
                right_order_mask = []
                for el in x_ordering:
                        right_order_mask.append(np.where(x == el)[0][0])
                
                x = x[right_order_mask]
                y = y[right_order_mask]
                y_err = y_err[right_order_mask]

        fig, ax = plt.subplots(figsize=plot_params['figsize'])
        ax.grid(axis = grid_axis, 
                linewidth = grid_linewidth,
                linestyle = grid_linestyle)

        ax.set_title(title, fontsize=plot_params['label_size'])

        ax.plot(x, y, markersize = plot_params['markersize'],
                marker='o', mec=color, 
                mfc='none', linestyle='None')
        ax.errorbar(x, y, yerr=y_err, lw=plot_params['linewidth'],
                fmt= '-', color = color)

        ax.tick_params(axis='both', which='major', 
                        labelsize=plot_params['label_size'])  

        ax.tick_params(axis = 'x', rotation=plot_params['xlabel_rotation'])  

        ax.set_xlabel(xlabel, fontsize = plot_params['label_size'])
        ax.set_ylabel('ID-score', 
                        fontsize = plot_params['label_size'])

        y_ticks = np.array([0, 0.25, 0.5, 0.75, 1])
        ax.set_yticks(y_ticks)
        y_tick_labels = [0.0, "", 0.5, "", 1.0]
        ax.set_yticklabels(y_tick_labels)

        sns.despine()




def plot_lineage_trend(ax, x_val, x, y, y_err, label, 
                       x_ordering, dataset_name, title = '',
                       plot_params = {"label_size": 25,
                                "markersize": 20,
                                "linestyle": ' ',
                                "linewidth": 4,
                                "alpha": 0.9, 
                                "markeredgewidth": 2,
                                "ticks_length": 2,
                                "ticks_width": 0.5,
                                "f_legend": 1}, 
                       text_cts = False):

        right_order_mask = []
        for el in x_ordering:
                right_order_mask.append(np.where(x == el)[0][0])
        
        x = x[right_order_mask]
        y = y[right_order_mask]
        y_err = y_err[right_order_mask]

        if(dataset_name != 'MousePancreas'):
                ax.set_title(title, fontsize=plot_params['label_size'])

        ax.plot(x_val, y, markersize = plot_params['markersize'],
                marker=plot_params['marker'], alpha = plot_params['alpha'], 
                mec=plot_params['color'], mfc='none', 
                mew = plot_params['markeredgewidth'],
                linestyle='None', label = label)
        
        if(plot_params['linestyle'] != ' '):
                ax.plot(x_val, y, '-', color = plot_params['color'], 
                        alpha = plot_params['alpha'], linewidth = plot_params['linewidth'], 
                        linestyle = plot_params['linestyle'])
        ax.errorbar(x_val, y, yerr=y_err, lw=plot_params['markeredgewidth'], 
                    alpha = plot_params['alpha'], fmt= ' ', 
                    color = plot_params['color'])

        ax.set_xlabel('', fontsize = plot_params['label_size'])
        ax.set_ylabel('ID-score', fontsize = plot_params['label_size'])
        ax.set_xlabel('<-----Cell potency----', fontsize = plot_params['label_size'])

        y_ticks = np.array([0, 0.25, 0.5, 0.75, 1])
        ax.set_yticks(y_ticks)
        y_tick_labels = [0.0, "", 0.5, "", 1.0]
        ax.set_yticklabels(y_tick_labels)

        ax.set_xticklabels([])

        legend_size = plot_params['label_size']*plot_params['f_legend']
        ax.legend(fontsize = legend_size, ncol = 1,
                fancybox=True, shadow = True,
                labelspacing = 0.5, loc = 'lower left',
                columnspacing=0.)
        
        ax.tick_params(axis = 'both', which = 'major', 
                        labelsize = plot_params['label_size'], 
                        length = plot_params['ticks_length'],
                        width = plot_params['ticks_width'])  

        f_shift = 0.04
        x_pos_ct = x_val*(1+f_shift)
        y_pos_ct = y*(1+f_shift)
        if('Endocrine' in x):   
                x_pos_ct[-1] = x_val[-1]*(1-4.7*f_shift)
                y_pos_ct[-1] = y[-1]*(1-15*f_shift)
        if('Tip' in x):   
                x_pos_ct[1] = x_val[1]*(1+3*f_shift)
                y_pos_ct[1] = y[1]*(1-0.75*f_shift)
        if('EP' in x):   y_pos_ct[2] = y[2]*(1-3.5*f_shift)
        if('Acinar' in x):   y_pos_ct[2] = y[2]*(1-1.5*f_shift)
        if('Fev+' in x):   y_pos_ct[3] = y[3]*(1-f_shift)
        if('Ductal' in x):   y_pos_ct[2] = y[2]*(1-0.65*f_shift)
        if('Multipotent' in x): x_pos_ct[0] = 2.5*f_shift
        
        text_size = plot_params['f_legend']*plot_params['label_size']
        if(text_cts):
                for n, ct_name in enumerate(x):
                        ax.text(x_pos_ct[n], y_pos_ct[n], ct_name, 
                                fontsize = text_size)
        
        if(dataset_name != 'MousePancreas'):
                ax.set_xticks([])
                ax.grid(axis = grid_axis, 
                        linewidth = grid_linewidth,
                        linestyle = grid_linestyle)
        
        sns.despine()




def tree_assembly(ax, x, y, y_err, dataset_name, 
                  title = '',
                  plot_params = {"label_size": 25,
                                "markersize": 20,
                                "linewidth": 4,
                                "alpha": 0.9, 
                                "markeredgewidth": 2,
                                "ticks_length": 2,
                                "ticks_width": 0.5,
                                "f_legend": 1}):

        lineage = ['Multipotent', 'Tip', 'Acinar']
        plot_params['color'] = 'deepskyblue'
        plot_params['linestyle'] = 'dashed'
        plot_params['marker'] = '^'
        label = 'Lineage 1'
        x_val = np.linspace(0, 10*1/2, len(lineage))
        plot_lineage_trend(ax, x_val, x, y, y_err, label, 
                           lineage, dataset_name,
                           plot_params=plot_params,
                           text_cts=True)

        lineage = ['Multipotent', 'Trunk', 'Ductal']
        plot_params['color'] = 'turquoise'
        plot_params['linestyle'] = 'solid'
        plot_params['marker'] = '>'
        label = 'Lineage 2'
        x_val = np.linspace(0, 10*1/2, len(lineage))
        plot_lineage_trend(ax, x_val, x, y, y_err, label, 
                           lineage, dataset_name,
                           plot_params=plot_params,
                           text_cts=True)

        lineage = ['Multipotent', 'Trunk', 'EP', 'Fev+', 'Endocrine']
        plot_params['color'] = 'deeppink'
        plot_params['linestyle'] = 'dotted'
        plot_params['marker'] = '<'
        label = 'Lineage 3'
        x_val = np.linspace(0, 10, len(lineage))
        plot_lineage_trend(ax, x_val, x, y, y_err, label, 
                           lineage, dataset_name,
                           plot_params=plot_params,
                           text_cts=True)
        
        ax.set_title('Mouse pancreatic endocrinogenesis\ndifferentiation tree', 
                     fontsize=plot_params['label_size'])
        
        ax.set_xticks([])

        ax.grid(axis = grid_axis, 
                linewidth = grid_linewidth,
                linestyle = grid_linestyle)

        sns.despine()


