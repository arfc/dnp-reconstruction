from scipy.integrate import cumtrapz
import numpy as np
import matplotlib.pyplot as plt
import os


def delnu_per_fiss(times, counts, fissions, efficiency):
    """
    Calculates the number of delayed neutrons per fission
        for a given dataset
    """
    norm_cnts = list()
    for cnt in counts:
        norm_cnts.append(cnt / (fissions * efficiency))
    #int_cnt = cumtrapz(norm_cnts,
    #                   x=times)
    #tot_cnt = int_cnt[-1] - int_cnt[0]
    tot_cnt = np.trapz(norm_cnts, x=times)
    dnpf = tot_cnt# / (fissions * efficiency)
    print(f'Approximate n/f: {dnpf}\n')
    print('-'*40)
    return dnpf

def delnu_per_fiss_norm(times, norm_counts):
    """
    Calculates the number of delayed neutrons per fission
        for a given dataset without eff or fiss
    """
    int_cnt = cumtrapz(norm_counts,
                       x=times)
    tot_cnt = int_cnt[-1] - int_cnt[0]
    dnpf = tot_cnt
    print(f'Approximate2 n/f: {dnpf}\n')
    print('-'*40)
    return dnpf

def multplt(x,
            y,
            label='unlabeled',
            alpha=1,
            errors=None):
    """
    Makes multiplotting on a single figure easier
    """
    linestyle_tuple = [
     ('dotted',                (0, (1, 1))),

     ('loosely dashed',        (0, (5, 10))),
     ('dashed',                (0, (5, 5))),

     ('loosely dashdotted',    (0, (3, 10, 1, 10))),
     ('dashdotted',            (0, (3, 5, 1, 5))),

     ('dashdotdotted',         (0, (3, 5, 1, 5, 1, 5))),
     ('loosely dashdotdotted', (0, (3, 10, 1, 10, 1, 10)))]
    #linestyle=linestyle_tuple[np.random.randint(0, len(linestyle_tuple))][-1]
    global miscmultpltid
    try:
        miscmultpltid
    except NameError:
        miscmultpltid = 0
    miscmultpltid += 1
    if miscmultpltid == len(linestyle_tuple):
        miscmultpltid = 0
    linestyle=linestyle_tuple[miscmultpltid][-1]
    plt.plot(x, y, label=label, alpha=alpha, linestyle=linestyle)
    if type(errors) != type(None):
        if type(errors) == type(list()):
            errors = np.array(errors)
        while len(np.shape(errors)) > 1:
            errors = errors[0]
        if type(y) == type(list()):
            y = np.array(y)
        while len(np.shape(y)) > 1:
            y = y[0]
        #plt.errorbar(x, y, yerr=errors, elinewidth=1,
                     #alpha=alpha, linestyle=linestyle, label=label)
        plt.fill_between(x, y+errors, y-errors, alpha=alpha/2)
        
    else:
        #
        pass
    plt.legend()
    return

def dir_handle(path):
    """
    Checks if the directory exists, and creates it if it doesn't
    """
    if not os.path.isdir(path):
        val = os.mkdir(f'{path}')
    return
    

def movie_gen(current_path, num_files):
    """
    Generate a .gif

    Parameters
    ----------
    current_path : str
        Current path where files are located
    num_files : int
        Number of plots

    Returns
    -------
    None
    """

    filenames = list()
    for i in range(num_files):
        name = f'{i}.png'
        filenames.append(current_path + name)

    import imageio
    images = []
    for filename in filenames:
        images.append(imageio.v2.imread(filename))
    imageio.mimsave(f'{current_path}movie.gif', images) 

    return
