from mosden.postprocessing import PostProcess
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

class MultiPostProcess():
    def __init__(self, input_paths: list[str]) -> None:
        self.posts: list[PostProcess] = [PostProcess(p) for p in input_paths]
        self.heatmap_key: str = 'modeling_options'
        self.heatmap_x: str = 'incore_s'
        self.heatmap_y: str = 'excore_s'
        self.hm_x_name = r'$\tau_{in}$ $[s]$ '
        self.hm_y_name = r'$\tau_{ex}$ $[s]$'
        self.hm_x_vals: list = list()
        self.hm_y_vals: list = list()
        self._initialize_posts()
        self.hm_z_names: dict[str, str] = {
            'summed_yield': r'$\bar{\nu}_d$',
            'group_yield': r'$\bar{\nu}_d$',
            'summed_avg_halflife': r'$\bar{T} [s]$',
            'group_avg_halflife': r'$\bar{T} [s]$'
        }
        return None
    
    def _initialize_posts(self):
        for post in self.posts:
            post.run()
            modeling_options: dict = post.input_data.get(self.heatmap_key, {})
            post.hm_x = modeling_options.get(self.heatmap_x, 0.0)
            post.hm_y = modeling_options.get(self.heatmap_y, 0.0)
            self.hm_x_vals.append(post.hm_x)
            self.hm_y_vals.append(post.hm_y)
        return None
    
    def run(self):
        self.heatmap_gen()
        self.group_param_histogram()
        return None
    
    def _collect_post_data(self) -> dict[str: list[float]]:
        post_data = dict()
        summed_yield = list()
        summed_avg_halflife = list()
        group_yield = list()
        group_avg_halflife = list()
        for post in self.posts:
            summed_yield.append(post.summed_yield.n)
            summed_avg_halflife.append(post.summed_avg_halflife.n)
            group_yield.append(post.group_yield.n)
            group_avg_halflife.append(post.group_avg_halflife.n)
        post_data['summed_yield'] = summed_yield
        post_data['summed_avg_halflife'] = summed_avg_halflife
        post_data['group_yield'] = group_yield
        post_data['group_avg_halflife'] = group_avg_halflife
        return post_data


    def heatmap_gen(self):
        z_values: dict[str, list[float]] = self._collect_post_data()
        for z_id in self.hm_z_names.keys():
            X = self.hm_x_vals
            Y = self.hm_y_vals
            z_name = self.hm_z_names[z_id]
            Z = z_values[z_id]

            df = pd.DataFrame.from_dict(np.array([X, Y, Z]).T)
            df.columns = [self.hm_x_name, self.hm_y_name, z_name]
            df[z_name] = pd.to_numeric(df[z_name])
            pivotted = df.pivot(index=self.hm_x_name, columns=self.hm_y_name, values=z_name)
            color = sns.color_palette("dark:pink_r", as_cmap=True)
            ax = sns.heatmap(pivotted, cmap=color)
            ax.invert_yaxis()
            ax.collections[0].colorbar.set_label(z_name)
            plt.tight_layout()
            plt.savefig(f'surf_{z_id}.png')
            plt.close()
        return None
    
    def group_param_histogram(self):
        data = dict()
        yields = list()
        halflives = list()
        groups = list()
        names = list()
        for post in self.posts:
            yield_val, halflife_val = post._get_MC_group_params()
            for gi, val in enumerate(yield_val):
                yields.append(val[0])
                halflives.append(halflife_val[gi][0])
                groups.append(gi+1)
                names.append(post.name)
        data['Yield'] = yields
        data['Halflife [s]'] = halflives
        data['Group'] = groups
        data['Data Source'] = names
        df = pd.DataFrame(data)
        df = df.explode(['Yield', 'Halflife [s]', 'Group'], ignore_index=True)
        sns.barplot(df, x='Group', y='Yield', hue='Data Source')
        plt.legend()
        plt.tight_layout()
        plt.savefig(f'yields.png')
        plt.close()
        sns.barplot(df, x='Group', y='Halflife [s]', hue='Data Source')
        plt.legend()
        plt.yscale('log')
        plt.tight_layout()
        plt.savefig(f'halflives.png')
        plt.close()
