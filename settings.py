import time
import numpy as np

begin = time.time()
fissions = 1.46516E18 # W17x17
fissions_error = 6.464567E13
volume = 108961.423 # W17x17
efficiency = 3.754492213194589e-08
mass_normalize = 1
alpha = 0.7
show_iso = 'all'
TRITON_out = './scale_outputs/godiva_3d_depl.out'
ensdf_fname = './ensdf_data/eval_net.xlsx'
ensdf_sheet = 'Sheet1'
endf_spectra_filename = './spectra/ENDF.xlsx'
endf_spectra_sheetname = 'BVII'


# Ensure this is corret
sample = 'uranium' # e50_dec_3 or e50_dec_33 or uranium or plutonium
irradiation = 'pulse' # pulse or infinite
input(f'Problem: {sample} with {irradiation} irradiation [Enter to Continue]')
print('-'*100)


# Common Options
dt = 1
end_time = 330 #2
start_time = 0
use_errorbars = True
target = 'all'
decay_nodes = 3
percent_variation = [0.1, 0.05, 0.04, 0.03, 0.02, 0.01, 0.005, 0.001]

DEBUG_IGNORE_ISOTOPES = []# ['ge86', 'i140', 'y98m', 'y97', 'i137', 'br93']

fit_groups = 6
decay_dims = 2
abund_iters = 100
spectra_iters = 100
#halflife_base = np.array([54.51, 21.84, 6.00, 2.23, 0.496, 0.179]) # Keepin
#halflife_base = np.array([55.64, 24.59, 4.338, 2.704, 0.378, 0.222]) # Important Isos
#halflife_base = np.array([53.14725, 21.27271, 5.83048, 2.19182, 0.52832, 0.11451]) # ORI-BEST
#halflife_base = np.array([52.12799, 22.54858, 4.74958, 1.83765, 0.45254, 0.09341]) # ORI-BEST Pu
halflife_base = np.array([48.79745, 19.25847, 3.6288, 1.26018, 0.31798, 0.09831]) # IAEA-BEST
#halflife_base = np.array([53.32415, 23.066, 6.36757, 2.38545, 0.95107, 0.30229]) # IAEA-BEST Pu
#halflife_base = np.array([53.6869, 21.51022, 4.91714, 2.01612, 0.6467, 0.13188]) # e50_dec_33 ORI-BEST
#halflife_base = np.array([58.65942, 25.76584, 11.83736, 3.44084, 11.83325, 5.28915]) # e50_dec_3 IAEA-BEST
#halflife_base = np.array([49.0, 19.2, 3.64, 1.28, 0.320, 0.098])
#halflife_base = np.array([51.3, 20.7, 6.04, 2.19, 0.505, 0.115])
times = np.arange(start_time, end_time+dt, dt)

# For good results, run with many fine time steps (go out to inf)
run_compare_decay   = False # ORIGEN and IAEA compare decay constants
run_ori_tri_compare = False # ORIGEN and TRITON compare
run_keep_brad_scale = False # Keepin, Brady England, SCALE
ori_pure_ensdf_comp = False # Compare Pure ORIGEN and IAEA ORIGEN
ori_ensdf_counts_cmp= False # Compare IAEA ORIGEN with Pure ORIGEN
# Misc runs
collect_data        = False # Print total yields for different approaches
ori_ensdf_keep_err  = False # Compare uncertainties for ORIGEN, IAEA, and Keepin
test_custom_fit     = False # Use custom group (in results_gen)
triton_no_ori_over  = False # Show isotopes in TRITON but not in ORIGEN
keepin_pure         = False # Keepin data
find_worst_lam_isos = False # Isotopes biggest diff IAEA/ORI lambdas
find_worst_pn_isos  = True # Isotopes biggest diff IAEA/ORI Pn
find_worst_tot_isos = False # Isotopes biggest diff IAEA/ORI net
view_pn_ori_pure    = False # ORIGEN with IAEA Pn values
iaea_ori_pure       = False # Compare IAEA ORIGEN and Pure ORIGEN
important_pure_ori  = False # Print the most impactful isotopes pure ORIGEN
keepin_pure_ori_fit = False # Compare Pure ORIGEN fit with Keepin/BradyEngland
sanity_check        = False  # Keepin group fit and Keepin data 
yield_contributions = False  # Print out isotopic yield contributions
# Requests (7/13/22)
keepin_data_origen  = False # Keepin counts vs ORIGEN counts
keepin_base_compare = False # Keepin normalized group fits with B/E and ORIGEN group fits
integral_keep_ori_be= False # Integral as a function of time for the different group fits

# Targeted
# Lambda
targets_iso_iaea_ori = False # Targets for IAEA with ORIGEN decays
#target_list_iaea_ori = ['ge86', 'rb95', 'br90', 'rb94', 'br89'] # Old
target_list_iaea_ori = ['ge86', 'y98m', 'i137', 'i140', 'y97'] # New
# Pn
targets_pure_ori_ori = False # Targets for IAEA with ORIGEN Pns
target_list_ori_ori  = ['ge86', 'br91', 'y100', 'as85', 'i138'] # Old
target_list_ori_ori  = ['ge86', 'br91', 'as86', 'as85', 'i138'] # New

# Net
targets_net_pure_ori = False # Targets for Pure ORIGEN
#target_list_pure_tot = ['br91', 'ge86', 'as85', 'i138', 'br88'] # Short lived
#target_list_pure_tot = ['br88', 'br87', 'i137', 'te136', 'cs141'] # Old
target_list_pure_tot = ['br91', 'as86', 'as85', 'i138', 'i137'] # New
#target_list = ['all']

# For good results, run with coarse time steps (330 final, 1 step) (more nodes better)
test_group_fit      = False # Use test group (see results_gen)
ori_ensdf_group_fit = False #True # IAEA ORIGEN group fit
tri_ensdf_group_fit = False # IAEA TRITON group fit
ori_pure_group_fit  = False # Pure ORIGEN group fit
group_abundance_err = True #True # Default on, calculate errors using stochastic method


# PRKEs
reactivity_magnitudes = [0.5]
gen_time = 1E-7 #5.56122E-9 #1E-5
prk_dt = 1E-5
prk_tf = 1
nubar = 2.60340
prk_times = np.arange(0, prk_tf+prk_dt, prk_dt)
ori_iaea_keepin_prk = False # IAEA ORIGEN with Keepin PRK response
puori_iaea_ori_prk  = False # Pure ORIGEN with IAEA ORIGEN response
keepin_pure_iaea_prk= False # Keepin, Pure, and IAEA ORIGEN response


# Spectras
energy_mesh = np.linspace(0, 1.8, 10) #500I 1.8e6 1e-5
spectra_normalized  = True # Normalizes by dividing by sum of energy bins at given time; probability (default True)
spectra_uncertainty = False # Plot with or without uncertainties (negligible)
pure_ori_t_spectra  = False # Pure ORIGEN spectral results over time
pure_ori_2d_spectra = False # Pure ORIGEN spectral results in 2D matrix form
iaea_ori_t_spectra  = False # IAEA ORIGEN spectral results over time
iaea_ori_2d_spectra = False # IAEA ORIGEN spectral results in 2D matrix form
spectra_puoriaea_com= False # Comparison of IAEA and Pure ORIGEN spectra
display_endf_spectra= False # Generate plots of ENDF group spectra

# Spectra Fitting
spectra_expstrp_oria= False # Generate group spectra using exponential stripping
spectra_lstsq_oriaea= False # IAEA ORIGEN group spectra using fraction fitting
alt_spc_lstsq_oriaea= False # IAEA ORIGEN group spectra using data least squares
spec_compare_oriaea = False # Compare fraction fitting and data least squares
spectra_puori_fit   = False # Pure ORIGEN group spectra using data least squares



if sample == 'uranium':
    print('Using uranium sample')
    ORIGEN_out = './scale_outputs/godiva_irrad_post_pulse.out'
    imdir = './images/'
    fissions = 1.013343827616795e+16 # Godiva Pulse U Sample
    volume = 0.1583105694 # Godiva Pulse Volume
    mass_normalize  = 21.90177 # Godiva Pulse Mnorm
    # Data
    # Pure ORIGEN 0.281% diff
    pure_ori_lamvec = np.log(2) / np.array([53.14725, 21.27271, 5.83048, 2.19182, 0.52832, 0.11451])
    pure_ori_lamerr = np.log(2) / (np.array([53.14725, 21.27271, 5.83048, 2.19182, 0.52832, 0.11451]))**2 * np.array([0.266, 0.106, 0.03, 0.01, 0.003, 0.0006])
    pure_ori_abuvec = np.array([0.064106, 0.311083, 0.263024, 0.676943, 0.309569, 0.092138]) / 100
    pure_ori_abuerr = np.array([1.9044815832793787e-05, 3.1672704807472236e-05, 4.199253210280293e-05, 3.865434248852521e-05, 1.5348959476351932e-05, 5.038437840935975e-06])

    # ORIGEN-IAEA 0.566% diff
    ori_iaea_lamvec = np.log(2) / np.array([48.79745, 19.25847, 3.6288, 1.26018, 0.31798, 0.09831])
    ori_iaea_lamerr = np.log(2) / (np.array([48.79745, 19.25847, 3.6288, 1.26018, 0.31798, 0.09831]))**2 * np.array([0.976, 0.385, 0.073, 0.03, 0.006, 0.002])
    ori_iaea_abuvec = np.array([0.00082825, 0.00347624, 0.00673826, 0.00560417, 0.00181652, 0.00051663])
    ori_iaea_abuerr = np.array([0.00038698359751450075, 0.000407837027550829, 0.0002061831165613984, 0.00014827854962492612, 4.101993906021232e-05, 2.048824039853155e-05])

    # Keepin
    keepin_lamvec = np.log(2) / np.array([54.51, 21.84, 6, 2.23, 0.496, 0.179])
    keepin_lamerr = (np.log(2) /
                              (np.array([54.51, 21.84, 6, 2.23, 0.496, 0.179]))**2 *
                              np.array([0.94, 0.54, 0.17, 0.06, 0.029, 0.017]))
    keepin_abuvec = np.array([0.00063, 0.00351, 0.00310, 0.00672, 0.00211, 0.00043])
    keepin_abuerr = np.array([0.00005, 0.00011, 0.00028, 0.00023, 0.00015, 0.00005])

    # Brady England
    be_lamvec = np.log(2) / np.array([52.116, 21.197, 5.7380, 2.2891, 0.8159, 0.2430])
    be_lamerr = np.log(2) / (np.array([52.116, 21.197, 5.7380, 2.2891, 0.8159, 0.2430]))**2 * np.array([0, 0, 0, 0, 0, 0])
    be_abuvec = np.array([0.000721, 0.00372242, 0.0035535, 0.00796808, 0.00326716, 0.00136784])
    be_abuerr = np.array([1.442E-6, 7.4E-6, 7.1E-6, 1.6E-5, 6.5E-6, 2.7E-6])


    # SCALE (scale-beta.h5) 
    # 1 fast; 2 thermal
    scale_yield = 0.0064 * nubar # Roughly 0.1666
    scale_lamvec = np.array([0.0127, 0.0317, 0.115, 0.311, 1.4, 3.87])
    scale_abuvec = np.array([0.038, 0.213, 0.188, 0.407, 0.128, 0.026]) * scale_yield


elif sample == 'plutonium':
    print('Using plutonium sample')
    ORIGEN_out = './scale_outputs/godiva_irrad_post_pulse_pu.out'
    imdir = './pu-images/'
    fissions = 1.3709773371924748e+16 # Godiva Pulse Pu Sample
    volume = 0.1583105694 # Godiva Pulse Volume
    mass_normalize  = 21.90177 # Godiva Pulse Mnorm
    # Data
    # Pure ORIGEN 0.370% diff
    pure_ori_lamvec = np.log(2) / np.array([52.12799, 22.54858, 4.74958, 1.83765, 0.45254, 0.09341])
    pure_ori_lamerr = np.log(2) / (np.array([52.12799, 22.54858, 4.74958, 1.83765, 0.45254, 0.09341]))**2 * np.array([0.2606399499999945, 0.11274289999999887, 0.0237478999999996, 0.009188249999999898, 0.0022626999999999786, 0.00046704999999999663])
    pure_ori_abuvec = np.array([2.62248342e-04, 1.90035919e-03, 1.54918642e-03, 2.16181983e-03, 8.72871090e-04, 4.34401108e-05])
    pure_ori_abuerr = np.array([1.2026090153137625e-05, 1.643066269592019e-05, 1.7086814964549477e-05, 1.610957561996513e-05, 5.442227200279379e-06, 5.250137670049115e-07])

    # ORIGEN-IAEA 0.153% diff
    ori_iaea_lamvec = np.log(2) / np.array([53.32415, 23.066, 6.36757, 2.38545, 0.95107, 0.30229])
    ori_iaea_lamerr = np.log(2) / (np.array([53.32415, 23.066, 6.36757, 2.38545, 0.95107, 0.30229]))**2 * np.array([0.05332414999999813, 0.023065999999998255, 0.006367569999999656, 0.002385449999999789, 0.0009510699999999983, 0.00030228999999998285])
    ori_iaea_abuvec = np.array([0.00024107, 0.00195117, 0.0012267,  0.00244866, 0.00075103, 0.0004652 ])
    ori_iaea_abuerr = np.array([7.488384862926156e-06, 2.3910588399809653e-05, 1.1354276614836517e-05, 3.486223251542677e-05, 7.426358566914596e-06, 9.901705477276028e-06])


    # Keepin
    keepin_lamvec = np.log(2) / np.array([53.75, 22.29, 5.19, 2.09, 0.549, 0.216])
    keepin_lamerr = (np.log(2) /
                              (np.array([53.75, 22.29, 5.19, 2.09, 0.549, 0.216]) )**2 *
                              np.array([0.95, 0.36, 0.12, 0.08, 0.049, 0.017]))
    keepin_abuvec = np.array([0.024, 0.176, 0.136, 0.207, 0.065, 0.022]) / 100
    keepin_abuerr = np.array([0.002, 0.009, 0.013, 0.012, 0.007, 0.003]) / 100

    # Brady England
    be_lamvec = np.array([0.0133, 0.0309, 0.1134, 0.2925, 0.8575, 2.7297])
    be_lamerr = np.array([0, 0, 0, 0, 0, 0])
    be_abuvec = np.array([0.0363, 0.2364, 0.1789, 0.3267, 0.1702, 0.0515]) * 0.68/100
    be_abuerr = be_abuvec * 0.08/100

else:
    print(f'Using {sample} sample')
    ORIGEN_out = f'./scale_outputs/{sample}.out'
    imdir = f'./{sample}-images/'


# IAEA Spectra File Names and associated isotope
# Using rb97 spectrum for rb98 (very similar)
# br92 only has figure, no data
iaea_spectra = {'b15' : ['B-15-2000Bu33-fig11.txt',
                         'B-15-2003Mi01-fig6.txt'],
                'b17' : ['B-17-1997YaZX-fig3.txt'],
                'be14' : ['Be-14-1997Ao04-fig3.txt'],
                'c16' : ['C-16-2001Gr06-fig3.txt'],
                'cs141' : ['Cs-141-1989BRZI-fig39-full.txt'],
                'cs142' : ['Cs-142-1989BRZI-fig40-full.txt'],
                'cs143' : ['Cs-143-1989BRZI-fig41-full.txt'],
                'cs144' : ['Cs-144-1989BRZI-fig42-full.txt'],
                'cs145' : ['Cs-145-1989BRZI-fig43-full.txt'],
                'cs146' : ['Cs-146-1989BRZI-fig44-full.txt'],
                'cs147' : ['Cs-147-1989BRZI-fig45-full.txt'],
                'ga79' : ['fig12_full.dat'],
                'ga80' : ['fig13_full.dat'],
                'ga81' : ['fig14_full.dat'],
                'as85' : ['fig15_full.dat'],
                'br87' : ['fig16_full.dat'],
                'br88' : ['fig17_full.dat'],
                'br89' : ['fig18_full.dat'],
                'br90' : ['fig19_full.dat'],
                'br91' : ['fig20_full.dat'],
                #'br92' : ['fig21_full.dat'],
                'rb92' : ['fig22_full.dat'],
                'rb93' : ['fig23_full.dat'],
                'rb94' : ['fig24_full.dat'],
                'rb95' : ['fig25_full.dat'],
                'rb96' : ['fig26_full.dat'],
                'rb97' : ['fig27_full.dat'],
                'rb98' : ['fig27_full.dat'],
                'in129' : ['In-129-1989BRZI-fig29-full.txt'],
                'in130' : ['In-130-1989BRZI-fig30-full.txt'],
                'sn134' : ['Sn-134-1989BRZI-fig31-full.txt'],
                'sb135' : ['Sb-135-1989BRZI-fig32-full.txt'],
                'te136' : ['Te-136-1989BRZI-fig33-full.txt'],
                'i137' : ['I-137-1989BRZI-fig34-full.txt'],
                'i138' : ['I-138-1989BRZI-fig35-full.txt'],
                'i139' : ['I-139-1989BRZI-fig36-full.txt'],
                'i140' : ['I-140-1989BRZI-fig37-full.txt'],
                'i141' : ['I-141-1989BRZI-fig38-full.txt'],
                'he8' : ['He-8-1981Bj03-fig4.txt'],
                'li9' : ['Li-9-1990Ny01-fig2.txt',
                         'Li-9-2015Hi02-fig7.txt'],
                'li11' : ['Li-11-1979Az03-fig1.txt',
                          'Li-11-1979Az03-fig1_inset.txt',
                          'Li-11-1997Ao04-fig1.txt,'
                          'Li-11-1997Mo35-fig1.txt',
                          'Li-11-2004Hi24-fig1.txt'],
                'be15' : [''],
                'n17' : ['N-17-2000Bu33-fig12.txt',
                         'N-17-2001Gr06-fig1.txt',
                         'N-17-2003Mi01-fig5.txt'],
                'n18' : ['N-18-2007Lo05-fig2.txt'],
                'n21' : ['N-21-2008Lo06-fig4a.txt',
                         'N-21-2009Li51-fig6b.txt'],
                'na27' : ['Na-27-1981ZiZW-fig1.txt'],
                'na28' : ['Na-28-1981ZiZW-fig3.txt']
                }
endf_spectra = ['ge86', 'y98m', 'y99', 'br93']                



normalize_value = mass_normalize * volume




print('Running...')
if __name__ == '__main__':
    print('Not recommended to run settings script. Try `results_gen.py`')
    raise Exception
