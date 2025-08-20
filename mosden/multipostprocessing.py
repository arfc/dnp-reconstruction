from mosden.postprocessing import PostProcess
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from mosden.countrate import CountRate
import os

class MultiPostProcess():
    def __init__(self, input_paths: list[str]) -> None:
        self.posts: list[PostProcess] = [PostProcess(p) for p in input_paths]
        self.output_dir = self.posts[0].output_dir
        if len(self.posts) > 1:
            self.output_dir = f'./{self.posts[0].multi_id}/images/'
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

        self.do_heatmap = False


        self.hm_x_vals: list = list()
        self.hm_y_vals: list = list()
        self._post_heatmap_setup()
        self._initialize_posts()
        self.hm_z_names: dict[str, str] = {
            'summed_yield': r'$\bar{\nu}_d$',
            'group_yield': r'$\bar{\nu}_d$',
            'summed_avg_halflife': r'$\bar{T} [s]$',
            'group_avg_halflife': r'$\bar{T} [s]$'
        }
        self._set_post_names()
        return None
    
    def _set_post_names(self):
        if np.all([post.multi_id == 'tintex' for post in self.posts]):
            for post in self.posts:
                post.name = f'({post.t_in}, {post.t_ex})'
        elif np.all([post.multi_id == 'chem_long' for post in self.posts]):
            self.posts[0].name = 'Full MSBR'
            self.posts[1].name = 'Partial MSBR'
        elif np.all([post.multi_id == 'chem_bool' for post in self.posts]):
            self.posts[0].name = 'Full MSBR'
            self.posts[1].name = 'No removal'
        return None
    
    def _post_heatmap_setup(self):
        if np.all([post.multi_id == 'tintex' for post in self.posts]):
            self.heatmap_key: str = 'modeling_options'
            self.heatmap_x: str = 'incore_s'
            self.heatmap_y: str = 'excore_s'
            self.hm_x_name = r'$\tau_{in}$'
            self.hm_x_units = r'$[s]$'
            self.hm_y_name = r'$\tau_{ex}$'
            self.hm_y_units = r'$[s]$'
            self.do_heatmap = True
        return None
 
    def _initialize_posts(self):
        for post in self.posts:
            post.run()
            try:
                modeling_options: dict = post.input_data.get(self.heatmap_key, {})
                post.hm_x = modeling_options.get(self.heatmap_x, 0.0)
                post.hm_y = modeling_options.get(self.heatmap_y, 0.0)
                self.hm_x_vals.append(post.hm_x)
                self.hm_y_vals.append(post.hm_y)
            except AttributeError:
                self.do_heatmap = False
        return None
    
    def run(self):
        if self.do_heatmap:
            self.heatmap_gen()
        self.group_param_histogram()
        self.group_fit_counts()
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
            x_name = self.hm_x_name + self.hm_x_units
            X = self.hm_x_vals
            y_name = self.hm_y_name + self.hm_y_units
            Y = self.hm_y_vals
            z_name = self.hm_z_names[z_id]
            Z = z_values[z_id]

            df = pd.DataFrame.from_dict(np.array([X, Y, Z]).T)
            df.columns = [x_name, y_name, z_name]
            df[z_name] = pd.to_numeric(df[z_name])
            pivotted = df.pivot(index=x_name, columns=y_name, values=z_name)
            color = sns.color_palette("dark:pink_r", as_cmap=True)
            ax = sns.heatmap(pivotted, cmap=color)
            ax.invert_yaxis()
            ax.collections[0].colorbar.set_label(z_name)
            plt.tight_layout()
            plt.savefig(f'{self.output_dir}surf_{z_id}.png')
            plt.close()

            x_to_y_ratio = np.asarray(X)/np.asarray(Y)
            sorted_indices = np.argsort(x_to_y_ratio)
            x_to_y_ratio = np.array(x_to_y_ratio)[sorted_indices]
            sorted_z = np.array(Z)[sorted_indices]

            plt.plot(x_to_y_ratio, sorted_z, marker='.', markersize=5)
            plt.xlabel(f'{self.hm_x_name}/{self.hm_y_name}')
            plt.xscale('log')
            plt.ylabel(z_name)
            plt.tight_layout()
            plt.savefig(f'{self.output_dir}ratio_{z_id}.png')
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
        plt.savefig(f'{self.output_dir}yields.png')
        plt.close()
        sns.barplot(df, x='Group', y='Halflife [s]', hue='Data Source')
        plt.legend()
        plt.yscale('log')
        plt.tight_layout()
        plt.savefig(f'{self.output_dir}halflives.png')
        plt.close()
        return None
    
    def group_fit_counts(self):
        for pi, post in enumerate(self.posts):
            times = post.decay_times
            countrate = CountRate(post.input_path)
            countrate.method = 'groupfit'
            group_counts = countrate.calculate_count_rate(write_data=False)
            plt.plot(
                times,
                group_counts['counts'],
                color=f'C{pi}',
                alpha=0.75,
                label=post.name,
                linestyle='--',
                zorder=3,
                marker=post.markers[pi%len(post.markers)],
                markersize=3,
                markevery=5)
            plt.fill_between(
                times,
                group_counts['counts'] -
                group_counts['sigma counts'],
                group_counts['counts'] +
                group_counts['sigma counts'],
                color=f'C{pi}',
                alpha=0.3,
                zorder=2,
                edgecolor='black')
        plt.xlabel('Time [s]')
        plt.ylabel(r'Count Rate $[n \cdot s^{-1}]$')
        plt.yscale('log')
        plt.legend()
        plt.tight_layout()
        plt.savefig(f'{self.output_dir}group_counts.png')
        plt.close()


