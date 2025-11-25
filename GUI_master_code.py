import os, re
import sys
from multiprocessing import Pool

import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import numpy as np
from astropy.io import fits
from IPython.display import Image
from scipy.optimize import curve_fit
import rebound
from datetime import datetime
from rich import print
from datetime import datetime, timedelta
from dateutil.relativedelta import relativedelta

import pmoired
from pmoired import tellcorr

from rich.console import Console
from rich.table import Table
from rich.panel import Panel
from rich import box
import datetime as dt


# ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
# ┃                  Control Panel                     ┃
# ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
console = Console()

# Neon palette
C_FN      = "#ff5cd2" # function names (neon purple)
C_PARAM   = "#00ffcc" # parameters (neon cyan)
C_DEF     = "#ffae42" # defaults (neon orange)
C_FRAME   = "#ffae42" #"#39ff14"   # frame border (neon green)
ERR_COLOR = "#ff4747" # Error color red

def sig(fn, params):
    parts = []
    for p in params:
        if "=" in p:
            name, val = p.split("=")
            parts.append(f"[{C_PARAM}]{name}[/]=[{C_DEF}]{val}[/]")
        else:
            parts.append(f"[{C_PARAM}]{p}[/]")
    return f"[{C_FN}]{fn}[/](" + ", ".join(parts) + ")"

# Build the table
tbl = Table.grid(padding=(0, 2))
tbl.add_column("Description", justify="right", style="bright_black")
tbl.add_column("Function Signature")

tbl.add_row("Wavelength (bary + systemic corr).", sig("WAVE", []))
tbl.add_row("Extract raw normalized flux.......", sig("RawFlux", []))
tbl.add_row("Pollux spectra....................", sig("PolluxSpectra", []))
tbl.add_row("Norm Flux also from Pollux........", sig("NormFlux", []))
tbl.add_row("Total visibility amplitude........", sig("VISAMP", []))
tbl.add_row("Amplitude uncertainty.............", sig("VISAMPERR", []))
tbl.add_row("Closure phase (triple product)....", sig("T3PHI", []))
tbl.add_row("Differential phase................", sig("VISPHI", []))
tbl.add_row("Differential phase error..........", sig("VISPHIERR", []))
tbl.add_row("Total phase.......................", sig("VISPHITOT", []))
tbl.add_row("Continuum fit (amplitude/phase)...", sig("DATA_CONT", []))
tbl.add_row("Continuum differential phase......", sig("VISPHICONT", []))
tbl.add_row("Continuum visibility amplitude....", sig("VISAMPCONT", []))
tbl.add_row("Pure-line visibility..............", sig("VISLINE", []))
tbl.add_row("Mask Brγ wavelength range.........", sig("mask_brg", []))
tbl.add_row("Line-phase computation............", sig("PHILINE", []))
tbl.add_row("Print barycentric velocity info...", sig("print_bary", []))
tbl.add_row("Scatter epoch data point..........", sig("scatter_epoch_point", []))
tbl.add_row("Fit Gaussian to spectrum..........", sig("FitGaussian", []))
tbl.add_row("Annotate plot element.............", sig("annotate_func", []))
tbl.add_row("Orbit simulation and visualization", sig("orbit", []))
tbl.add_row("Theoretical differential phase....", sig("PHI_THEO", []))
tbl.add_row("Check for k values................", sig("apply_k", []))
tbl.add_row("Projected photocenter.............", sig("projected_p", []))
tbl.add_row("Fit photocenter from phase........", sig("p_phi_fit", []))

panel = Panel(
    tbl,
    title="[bold #ffae42]Available Functions[/]",
    border_style=C_FRAME,
    box=box.ROUNDED,
    padding=(1, 2)
)

#console.print(panel)



# ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
# ┃                      constants                     ┃
# ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
bg = 2.16612e-6 # meter
mas2rad = 4.848136811e-9
rad2mas = 1/mas2rad
c_km_s = 299792.458

# ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
# ┃                   Epoch number                     ┃
# ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

# in python convention with 0
file_num = {
    "2018_red":  7,
    "2019_red":  3,
    "2020_red":  3,
    "2021_red":  10,
    "2022a_red": 3,
    "2022b_red": 3,
    "2023_red":  4
}

# ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
# ┃                   Epoch number                     ┃
# ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
year_num = {
    "2018_red":  0,
    "2019_red":  1,
    "2020_red":  2,
    "2021_red":  3,
    "2022a_red": 4,
    "2022b_red": 5,
    "2023_red":  6
}

# ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
# ┃                    Flux ratios                     ┃
# ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
flux_ratios = {
    "2018_red":  0.520,
    "2019_red":  1.000,
    "2020_red":  0.429,
    "2021_red":  0.411,
    "2022a_red": 0.452,
    "2022b_red": 0.464,
    "2023_red":  0.367
}


# ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
# ┃           0 bracket gamma velocities               ┃
# ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

# --- Since radail velocity is difficult to estimate,
# --- We can zero the bracket gamma emission so it lies
# --- at the theoretical bracket gamma value.
# --- the variable is structured as:
# --- zeroBr[year][filenum] in km/s
# --- i.e. this is juts vmid
zeroBrg = {
    "2018_red":  [np.nan, 109.76, 49.95, 57.45, 53.79, 45.63, 67.66, 68.22],
    "2019_red":  [95.95, 93.38, 104.25, 97.16],
    "2020_red":  [44.81, 46.18, 49.19, 50.13],
    "2021_red":  [80.86, 87.57, 92.16, 91.83, 91.24, 94.29, 94.08, 90.27, 89.41, 91.02, 89.01],
    "2022a_red": [90.90, 92.65, 91.06, 88.45],
    "2022b_red": [94.79, 93.27, 95.06, 95.86],
    "2023_red":  [39.45, 47.4, 40.07, 46.33]
}

# ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
# ┃          Radial Velocities per epoch               ┃
# ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

radvel_primary = {
    "2018_red":  -14.356938,
    "2019_red":  15.150509,
    "2020_red":  8.0810061,
    "2021_red":  7.8621486,
    "2022a_red": 6.5666435,
    "2022b_red": 6.546166,
    "2023_red":  5.6713122
}


# ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
# ┃                  Astrometrics [mas]                ┃
# ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
x_astrometrics = np.array([-9.730,  0.011, -5.813, -6.218, -8.990, -9.051, -11.256])
y_astrometrics = np.array([-2.459, -0.017, -5.436, -5.602, -6.844, -6.862,  -7.797])

# corresponding 1σ uncertainties (mas)
x_err_astrometrics = np.array([0.05, 1.16, 0.03, 0.02, 0.06, 0.02, 0.01])
y_err_astrometrics = np.array([0.05, 1.80, 0.06, 0.02, 0.08, 0.02, 0.03])

# ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
# ┃                     Plot style                     ┃
# ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
plt.rcParams.update({
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "xtick.minor.visible": True,
    "ytick.minor.visible": True,
    "xtick.major.size": 4,
    "ytick.major.size": 4,''
    "legend.frameon": True,
    "grid.alpha": 0.25,
    "grid.linestyle": ":",
    "savefig.dpi": 600,
    "axes.grid": True,     # enables light grid like your upper plot
    "figure.dpi": 120,     # for notebook display
})


# ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
# ┃              Baseline name function                ┃
# ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

baseline_cache = {}  # (year, filenum) -> np.ndarray of labels

def _compute_baseline_labels(year, filenum):
    hdul = file_cache[year][filenum]
    arrays = {
        h.header['ARRNAME']: h.data
        for h in hdul
        if h.header.get('EXTNAME') == 'OI_ARRAY'
    }

    all_baselines = []

    for h in hdul:
        if h.header.get('EXTNAME') not in ('OI_VIS', 'OI_VIS2'):
            continue

        arrname   = h.header['ARRNAME']
        arr       = arrays[arrname]
        sta_names = arr['STA_NAME']
        sta_ids   = arr['STA_INDEX']
        sta_pairs = h.data['STA_INDEX']

        id_to_name = dict(zip(sta_ids, sta_names))
        bl = np.array([f"{id_to_name[i]}-{id_to_name[j]}" for (i, j) in sta_pairs])
        all_baselines.append(bl)

    unique_baselines = np.unique(np.concatenate(all_baselines))
    return unique_baselines[::-1]


def baseline_name(year, baseline):

    # baseline_name(year, filenum, baseline) is also possible,
    # but since the first file of each epoch has the same
    # configuration as all other files we can simplify with filenum = 0
    filenum = 0
    key = (year, filenum)
    if key not in baseline_cache:
        baseline_cache[key] = _compute_baseline_labels(year, filenum)
    labels = baseline_cache[key]
    return labels[baseline]


# ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
# ┃                 Baseline colors                    ┃
# ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
def baseline_colors(baseline):

    # --- some pretty colors uwu
    colors = [
        "#5ac4ee", "#63e3bf", "#b95cf4",
        "#8cf55a", "#f559b7", "#f7b54c"
    ]

    return colors[baseline]


# ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
# ┃                   Barycentrics                     ┃
# ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
# --- Barycentric correction
rv_2018 = np.array([-20.887, -20.907, -20.917, -20.935, -20.948, -20.965, -20.977, -20.993])
rv_2019 = np.array([-23.105, -23.127, -23.138, -23.158])
rv_2020 = np.array([7.256, 7.232, 7.193, 7.181])
rv_2021 = np.array([-13.099, -13.127, -13.141, -13.184, -13.211, -13.224, -13.252, -13.265, -13.291, -13.304, -13.328])
rv_2022a = np.array([-12.350, -12.377, -12.391, -12.406])
rv_2022b = np.array([-14.980, -15.006, -15.019, -15.032])
rv_2023 = np.array([-8.434, -8.463, -8.492, -8.506])

rv_epoch = {
    "2018":   rv_2018,
    "2019":   rv_2019,
    "2020":   rv_2020,
    "2021":   rv_2021,
    "2022a": rv_2022a,
    "2022b": rv_2022b,
    "2023":   rv_2023
}


# ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
# ┃                   Fast loading                     ┃
# ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
# Load the four FITS files once
file_cache = {}

def preload_year(year):
    d = f"All_Epochs/{year}/"
    files = sorted([os.path.join(d, f) for f in os.listdir(d) if f.endswith('.fits')])
    file_cache[year] = [fits.open(f, memmap=True) for f in files]

preload_year("2018_red")
preload_year("2019_red")
preload_year("2020_red")
preload_year("2021_red")
preload_year("2022a_red")
preload_year("2022b_red")
preload_year("2023_red")


# ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
# ┃                  Corrected Wavelength              ┃
# ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
def WAVE(
    year,
    filenum,
    SystemicVel=25.4,
    primaryRadvel=None,
    allignPollux=False
):

    # --- SystemicVel according to Beskrovnaya is 25.4 km/s instead of the 21.7 km/s from Simbad

    # --- Choose which file and epoch
    hdul = file_cache[year][filenum]

    # --- extract wavelength
    CORR_WAVE = hdul[13].data["CORR_WAVE"]

    if primaryRadvel is None:
        radvel_prim = radvel_primary[year]
    else:
        radvel_prim = primaryRadvel

    # --- vmid used to shift the pollux wavelength so it allign with the FLC peak
    if allignPollux:
        allignPollux = zeroBrg[year][filenum]
    else:
        allignPollux = 0

    # --- proper corretion: lam_corr =lam_data * (1 - vtot/c)
    vtot = rv_epoch[year[:-4]][filenum] + SystemicVel + radvel_prim #+ allignPollux
    return CORR_WAVE * (1 - vtot/c_km_s)



# ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
# ┃             Extract Data Fits file                 ┃
# ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
def Extract(year, filenum, dataType):

    hdul = file_cache[year][filenum]

    if dataType == "TELL_TRANS":
        return hdul[13].data["TELL_TRANS"]
    elif dataType == "CORR_SPEC":
        return hdul[13].data["CORR_SPEC"]
    elif dataType == "RAW_SPEC":
        return hdul[13].data["RAW_SPEC"]
    elif dataType == "CORR_CONT":
        return hdul[13].data["CORR_CONT"]
    elif dataType == "VISAMP":
        return hdul[9].data["VISAMP"]
    elif dataType == "VISAMPERR":
        return hdul[9].data["VISAMPERR"]
    elif dataType == "VISPHI":
        return hdul[9].data["VISPHI"]
    elif dataType == "VISPHIERR":
        return hdul[9].data["VISPHIERR"]
    elif dataType == "T3PHI":
        return hdul[11].data["T3PHI"]
    elif dataType == "UCOORD":
        return hdul[10].data["UCOORD"]
    elif dataType == "VCOORD":
        return hdul[10].data["VCOORD"]
    # --- Flux error is taken as the mean error on the
    # --- 4 telescopes. For simplicity.
    elif dataType == "FLUXERR":
        cont = hdul[13].data["CORR_CONT"]
        raw_error = np.mean(hdul[12].data["FLUXERR"], axis=0)
        cont      = hdul[13].data["CORR_CONT"]
        norm_err  = raw_error/cont
        return norm_err


# ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
# ┃                     Tellurics                      ┃
# ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
def Tellurics(year, filenum):
    return Extract(year=year, filenum=filenum, dataType="TELL_TRANS")


# ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
# ┃              Pollux spectra loading                ┃
# ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
def PolluxSpectra(
    year,
    filenum,
    showHelp=True,
    whichStar=None,
    whichSpectrum=None,
    fullResSpectra=True,
    plot=True,
    listSpectra=True,
    returnSpectra=True,
    returnHighResWave=False,
    returnLowResWave=False,
    includeParam=True

):

    # --- Show help
    if showHelp:

        tbl = Table.grid(padding=(0, 2))
        tbl.add_column("Description", justify="right", style="bright_black")
        tbl.add_column("Usage", style=C_PARAM)
        tbl.add_row("Disable this help........................", f"[{C_PARAM}]showHelp[/]=[{C_DEF}]False[/]")
        tbl.add_row("select year..............................", f"[{C_PARAM}]year[/]=[{C_DEF}]value[/]")
        tbl.add_row("select filenumber........................", f"[{C_PARAM}]filenum[/]=[{C_DEF}]value[/]")
        tbl.add_row("Choose Primary or secondary star spectrum", f"[{C_PARAM}]whichStar[/]=[{C_DEF}]”Primary/Secondary”[/]")
        tbl.add_row("list the available spextra...............", f"[{C_PARAM}]listSpectra[/]=[{C_DEF}]True[/]")
        tbl.add_row("Choose One of the available spectra......", f"[{C_PARAM}]whichSpectrum[/]=[{C_DEF}]value[/]")
        tbl.add_row("Use full resolution spectra..............", f"[{C_PARAM}]fullResSpectra[/]=[{C_DEF}]True[/]")
        tbl.add_row("Plot spectra (deafult=True)..............", f"[{C_PARAM}]plot[/]=[{C_DEF}]True[/]")
        tbl.add_row("return value for spectra.................", f"[{C_PARAM}]returnSpectra[/]=[{C_DEF}]True[/]")
        tbl.add_row("return High res wavelength...............", f"[{C_PARAM}]returnHighResWave[/]=[{C_DEF}]True[/]")
        tbl.add_row("return Low res wavelength (same WAVE())..", f"[{C_PARAM}]returnLowResWave[/]=[{C_DEF}]True[/]")
        tbl.add_row("include parameter box in plot............", f"[{C_PARAM}]includeParam[/]=[{C_DEF}]True[/]")

        panel = Panel(
            tbl,
            title="[bold #ffae42]Pollux Spectra Help[/]",
            border_style=C_FRAME,
            box=box.ROUNDED,
            padding=(1, 2)
        )

        console.print(panel)

    # --- make sure user is selecting a spectrum
    if whichStar is None:
        return console.print(f"select between primary or secondary star with [{C_PARAM}]whichStar[/]=[{C_DEF}]”Primary/Secondary”[/]")
    if whichSpectrum is None:
        return console.print(f"select which spectrum to use with [{C_PARAM}]whichSpectrum[/]=[{C_DEF}]number[/]")

    # --- Load Spectra
    # --- Extract data from fits data files
    FLUXNORM = Extract(year=year, filenum=filenum, dataType="CORR_SPEC")

    # --- Extract data from Pollux file
    file_cachePollux = {}
    dire = f"POLLUX/{whichStar + " spectra/"}"
    files = sorted([os.path.join(dire, f) for f in os.listdir(dire) if f.endswith('.FITS')])
    file_cachePollux["Primary.FITS"] = [fits.open(f, memmap=True) for f in files]
    hdulPollux = file_cachePollux["Primary.FITS"][whichSpectrum] # -> whichSpectrum should be a number


    # --- Since pollux is zeroed at the theoretical Brg emission we need to
    # --- shift it so it allign with the peak of the emission
    # --- proper corretion: lam_pollux_corr = lam_pollux_data * (1 - vtot/c)
    allignPollux = zeroBrg[year][filenum]
    vtot = - allignPollux
    WAVE_POLLUX     = hdulPollux[1].data["wavelength"]*1e-10 * (1 - vtot/c_km_s)
    WAVE_REF        = WAVE(year, filenum)
    FLUXNORM_POLLUX = hdulPollux[1].data["normalized flux"] # --> high res data


    # --- Reduce resolution to fit that of the data
    # --- FLUXNORM_POLLUX acts as alpha_i
    FLUXNORM_POLLUX_lowres = np.interp(WAVE_REF, WAVE_POLLUX, FLUXNORM_POLLUX)

    if returnHighResWave:
        return WAVE_POLLUX
    if returnLowResWave:
        return WAVE_REF

    if returnSpectra:

        # --- Choose if high res or low res is returned
        if fullResSpectra:
            return FLUXNORM_POLLUX
        else:
            return FLUXNORM_POLLUX_lowres

    # --- plots for the spectra
    if plot:

        # --- add box with pollux parameters
        if includeParam:
            # --- extract information of POLLUX model
            info_dir  = f"POLLUX/info {whichStar}"
            info_file = f"spectrum{whichSpectrum}.txt"

            info_path = os.path.join(info_dir, info_file)

            with open(info_path, "r", encoding="utf-8") as f:
                text = f.read()

            teff_match = re.search(r'Teff\s*=\s*[\'"]?(\d+)', text) # ---------- T_eff
            logg_match = re.search(r'logg\s*=\s*[\'"]?([\d.]+)', text) # ------- surface gravity
            turbvel_match = re.search(r'turbvel\s*=\s*[\'"]?([\d.]+)', text) # - turbulent velocity
            Mdot_match  = re.search(r"Mdot\s*=\s*'(-?[\d.]+)'", text) # -------- Mass loss
            Vinfty_match = re.search(r"Vinfty\s*=\s*'(-?[\d.]+)'", text) # ----- terminal velocity in km/s
            beta_match  = re.search(r"beta\s*=\s*'(-?[\d.]+)'", text) # -------- beta parameter

            # the velocity is given as:
            # v(r) = v_infty ( 1 - R_star/r(**beta))

            Teff = int(teff_match.group(1)) if teff_match else None
            logg = float(logg_match.group(1)) if logg_match else None
            turbvel = float(turbvel_match.group(1)) if turbvel_match else None
            Mdot   = float(Mdot_match.group(1)) if Mdot_match   else None
            Vinfty = float(Vinfty_match.group(1)) if Vinfty_match else None
            beta   = float(beta_match.group(1)) if beta_match   else None

            # format text block
            param_text = (
                f"$T_{{eff}}$   = {Teff:.0f} K\n"
                f"$\\log g$ = {logg:.2f}\n"
                f"$v_{{turb}}$  = {turbvel:.1f} km/s\n"
                f"$\\dot{{M}}$     = 10$^{{{Mdot:.2f}}}$ M$_\\odot$/yr\n"
                f"$v_\\infty$    = {Vinfty:.0f} km/s\n"
                f"$\\beta$      = {beta:.2f}"
            )

        # --- setup plot
        plt.figure(figsize=(6.3, 3.4))
        plt.xlim(2.16e-6, 2.172e-6)
        plt.ylim(0.75, 1.05)

        # --- high res spectrum
        if fullResSpectra:

            plt.title(f"Pollux deafult absorption spectrum {whichSpectrum}")
            plt.plot(
                WAVE_POLLUX,
                FLUXNORM_POLLUX,
                label="Synthetic spectrum"
            )

        # --- low res spectrum
        else:

            plt.title(f"Pollux adapted resolution spectrum {whichSpectrum} ({year}, File {filenum})")
            plt.plot(
                WAVE_REF,
                FLUXNORM_POLLUX_lowres,
                label="Synthetic spectrum"
            )
        if includeParam:
            plt.text(
                2.1602e-6, .89,  # position in data coordinates
                param_text,
                fontsize=9,
                va="top",
                ha="left",
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.8)
            )

        plt.xlabel("Wavelength [m]")
        plt.ylabel("Normalized flux")
        plt.legend(frameon=False, loc="upper right")
        plt.tight_layout()



# ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
# ┃                  Normalized Flux                   ┃
# ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
def NormFlux(
    year=None,
    filenum=None,
    showHelp=False,
    returnFLCcorr=False,
    returnFLC=False,
    returnFLCerr=False,
    PhotosphericCorr=False,
    alpha1_whichSpectrum=None,
    alpha2_whichSpectrum=None,
    alpha1=None,
    alpha2=None,
    beta1=None,
    beta2=None,
    returnPhotoCorrectionFlux=True,
    returnPollux=False,
    returnTellurics=True,
    returnBrg=True,
    returnRAWFLUX=True,
    returnAbsorbtion1=False,
    returnAbsorbtion2=False
):

    # --- return the normalized spectrum from PMOIRED
    if returnFLC and returnFLCcorr:
        return console.print("ERROR: Choose between returnFLC or returnFLCcorr not both")

    if returnFLC:
        return Extract(year=year, filenum=filenum, dataType="CORR_SPEC")

    if returnFLCerr:
        return Extract(year=year, filenum=filenum, dataType="FLUXERR")

# --- check how to apply alphas
    if PhotosphericCorr:

        # --- if both alphas arent None then fall back to default (1,1)
        if alpha1 is not None and alpha2 is not None:
            alpha1 = WAVE(year, filenum)/WAVE(year, filenum)
            alpha2 = WAVE(year, filenum)/WAVE(year, filenum)

        # --- only alpha1 set -> get Pollux for alpha2
        elif alpha1 is not None and alpha2 is None:
            alpha2 = PolluxSpectra(
                        year=year,
                        filenum=filenum,
                        whichStar="Secondary",
                        whichSpectrum=alpha2_whichSpectrum,
                        fullResSpectra=False,
                        showHelp=False
                        )

        # --- only alpha2 set -> get Pollux for alpha1
        elif alpha2 is not None and alpha1 is None:
            alpha1 = PolluxSpectra(
                        year=year,
                        filenum=filenum,
                        whichStar="Primary",
                        whichSpectrum=alpha1_whichSpectrum,
                        fullResSpectra=False,
                        showHelp=False
                        )

        # --- both None -> Pollux for both
        else:
            alpha1 = PolluxSpectra(
                        year=year,
                        filenum=filenum,
                        whichStar="Primary",
                        whichSpectrum=alpha1_whichSpectrum,
                        fullResSpectra=False,
                        showHelp=False
                        )

            alpha2 = PolluxSpectra(
                        year=year,
                        filenum=filenum,
                        whichStar="Secondary",
                        whichSpectrum=alpha2_whichSpectrum,
                        fullResSpectra=False,
                        showHelp=False
                        )

    # --- if PhotosphericCorr is False then default alphas (1,1)
    else:
        alpha1 = WAVE(year, filenum)/WAVE(year, filenum)
        alpha2 = WAVE(year, filenum)/WAVE(year, filenum)
        beta1  = 0
        beta2  = 0

    # betas default
    if beta1 is None:
        beta1 = 0
    if beta2 is None:
        beta2 = 0

    if returnAbsorbtion1:
        return alpha1

    if returnAbsorbtion2:
        return alpha2



##########################
    # --- Extract data from fits data files
    FLUXNORM = Extract(year=year, filenum=filenum, dataType="CORR_SPEC")

    WL            = WAVE(year, filenum)
    TELLURICS     = Tellurics(year, filenum)
    RAW_SPEC      = Extract(year, filenum, dataType="RAW_SPEC")
    CORR_CONT     = Extract(year, filenum, dataType="CORR_CONT")
    RAW_NOWM_FLUX = RAW_SPEC/CORR_CONT

    fc        = flux_ratios[year]
    term1     = (alpha1 + beta1) * (1 + beta2) + fc * (alpha2 + beta2) * (1 + beta1)
    term2     = (1 + fc) * (1 + beta1) * (1 + beta2)
    F_LC_CORR = FLUXNORM + 1 - term1/term2

    if returnFLCcorr:
        return F_LC_CORR


    if returnTellurics:
      return TELLURICS

    if returnRAWFLUX:
        return RAW_NOWM_FLUX

    if returnPhotoCorrectionFlux:
        return F_LC_CORR

    if returnAbsorbtion1:
        return alpha1

    if returnAbsorbtion2:
        return alpha1


# ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
# ┃                      Raw Flux                      ┃
# ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
def RawFlux(year, filenum):
    return Extract(year=year, filenum=filenum, dataType="RAW_SPEC")


# ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
# ┃                Visibilitiy AMP                     ┃
# ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
def VISAMP(year, filenum, baseline):
    return Extract(year=year, filenum=filenum, dataType="VISAMP")[baseline]

# ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
# ┃              Visibilitiy AMP ERR                   ┃
# ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
def VISAMPERR(year, filenum, baseline):
    return Extract(year=year, filenum=filenum, dataType="VISAMPERR")[baseline]


# ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
# ┃                   Closure Phase                    ┃
# ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
def T3PHI(year, filenum, tel_config):
    return Extract(year=year, filenum=filenum, dataType="T3PHI")[tel_config]


# ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
# ┃              Differential phase Raw                ┃
# ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
def VISPHI(year, filenum, baseline):
    return Extract(year=year, filenum=filenum, dataType="VISPHI")[baseline]


# ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
# ┃            Differential phase ERROR                ┃
# ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
def VISPHIERR(year, filenum, baseline):
    return Extract(year=year, filenum=filenum, dataType="VISPHIERR")[baseline]


# ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
# ┃             Total Differential Phase               ┃
# ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
def VISPHITOT(year, filenum, baseline):
    return VISPHI(year, filenum, baseline) - VISPHICONT(year, filenum, baseline)


# ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
# ┃                CONTINUUM MODEL                     ┃
# ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
def DATA_CONT(year, filenum, baseline, dataType):

    # --- Data extraction ---
    WL          = WAVE(year, filenum)
    if dataType == "VISAMP":
        DATA    = VISAMP(year, filenum, baseline)
        DATAERR = 0.02*VISAMP(year, filenum, baseline)
    elif dataType == "VISPHI":
        DATA = VISPHI(year, filenum, baseline)
        DATAERR = np.array([1 for i in range(len(VISPHI(year, filenum, baseline)))])
    else:
        return print("ERROR: Only dataType = \"VISAMP\" or \"VISPHI\" supported")

    # --- Only fit around Brγ and exclude the line channels ---
    lower_edge, upper_edge = 2.13e-6, 2.20e-6
    inner_remove_low, inner_remove_high = 2.161e-6, 2.172e-6

    valid = np.isfinite(DATA) & np.isfinite(DATAERR)

    mask = (
        (WL > lower_edge) & (WL < upper_edge)
        & ~((WL >= inner_remove_low) & (WL <= inner_remove_high))
        & valid
    )

    # --- weighted linear fit y = a*λ + b
    x = WL[mask]
    y = DATA[mask]
    w = 1.0 / np.clip(DATAERR[mask], 1e-12, np.inf)

    # coefficients + covariance
    (a, b), cov = np.polyfit(x, y, 1, w=w, cov=True)
    sa  = np.sqrt(cov[0, 0])
    sb  = np.sqrt(cov[1, 1])
    sab = cov[0, 1]

    # --- Define the fitted cont.
    xline = WAVE(year, filenum)
    cont_model = a * xline + b

    # --- residuals for future implementation
    yhat = a * x + b
    sigma = DATAERR[mask]
    residuals = (y - yhat) / sigma

    return cont_model


# ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
# ┃            Differential phase continuum            ┃ ✅
# ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
def VISPHICONT(year, filenum, baseline):
    return DATA_CONT(year, filenum, baseline, dataType="VISPHI")


# ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
# ┃                 VISAMP continuum                   ┃ ✅
# ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
def VISAMPCONT(year, filenum, baseline):
    return DATA_CONT(year, filenum, baseline, dataType="VISAMP")


# ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
# ┃                Pure line visibility                ┃ ✅
# ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
def VISLINE(
    year=None,
    filenum=None,
    baseline=None,
    PhotosphericCorr=False,
    alpha1_whichSpectrum=None,
    alpha2_whichSpectrum=None,
    alpha1=None,
    alpha2=None,
    beta1 =None,
    beta2 =None,
    showHelp=False,
    returnError=False
):

    # --- Show help
    if showHelp:

        tbl = Table.grid(padding=(0, 2))
        tbl.add_column("Description", justify="right", style="bright_black")
        tbl.add_column("Usage", style=C_PARAM)

        tbl.add_row("Disable this help..................", f"[{C_PARAM}]showHelp[/]=[{C_DEF}]False[/]")
        tbl.add_row("apply Photospheric correction......", f"[{C_PARAM}]PhotosphericCorr[/]=[{C_DEF}]True[/]")
        tbl.add_row("Return Error instead...............", f"[{C_PARAM}]returnError[/]=[{C_DEF}]True[/]")
        tbl.add_row("define alpha1......................", f"[{C_PARAM}]alpha1[/]=[{C_DEF}]True[/]")
        tbl.add_row("define alpha2......................", f"[{C_PARAM}]alpha2[/]=[{C_DEF}]True[/]")
        tbl.add_row("choose which alpha1 spectrum to use", f"[{C_PARAM}]alpha1_whichSpectrum[/]=[{C_DEF}]value[/]")
        tbl.add_row("choose which alpha2 spectrum to use", f"[{C_PARAM}]alpha2_whichSpectrum[/]=[{C_DEF}]value[/]")
        tbl.add_row("define beta1.......................", f"[{C_PARAM}]beta1[/]=[{C_DEF}]True[/]")
        tbl.add_row("define beta2.......................", f"[{C_PARAM}]beta2[/]=[{C_DEF}]True[/]")

        panel = Panel(
            tbl,
            title="[bold #ffae42]Pure Line visibility help[/]",
            border_style=C_FRAME,
            box=box.ROUNDED,
            padding=(1, 2)
        )

        console.print(panel)

    # --- Check if user inputed all necessary parameters
    if year is None or filenum is None or baseline is None:

        tbl = Table.grid(padding=(0, 2))
        tbl.add_column("Description", justify="right", style="bright_black")
        tbl.add_column("Usage", style=C_PARAM)

        tbl.add_row("choose Epoch.......", f"[{C_PARAM}]year[/]=[{C_DEF}]epoch[/]")
        tbl.add_row("choose File number.", f"[{C_PARAM}]filenum[/]=[{C_DEF}]value[/]")
        tbl.add_row("choose baseline....", f"[{C_PARAM}]baseline[/]=[{C_DEF}]value[/]")

        if PhotosphericCorr and (alpha1 is None or alpha2 is None):

            tbl.add_row("Specify alpha1.....................", f"[{C_PARAM}]alpha1[/]=[{C_DEF}]placeholder[/]")
            tbl.add_row("Specify alpha2.....................", f"[{C_PARAM}]alpha2[/]=[{C_DEF}]placeholder[/]")
            tbl.add_row("choose which alpha1 spectrum to use", f"[{C_PARAM}]alpha1_whichSpectrum[/]=[{C_DEF}]value[/]")
            tbl.add_row("choose which alpha2 spectrum to use", f"[{C_PARAM}]alpha2_whichSpectrum[/]=[{C_DEF}]value[/]")
            tbl.add_row("Specify beta1......................", f"[{C_PARAM}]beta1[/]=[{C_DEF}]placeholder[/]")
            tbl.add_row("Specify beta2......................", f"[{C_PARAM}]beta2[/]=[{C_DEF}]placeholder[/]")

        tbl.add_row("show help..........", f"[{C_PARAM}]showHelp[/]=[{C_DEF}]True[/]")

        panel = Panel(
            tbl,
            title="[bold #ff4747]ERROR --- Specify following parameters[/]",
            border_style=ERR_COLOR,
            box=box.ROUNDED,
            padding=(1, 2)
        )

        return console.print(panel)

    # --- Check if user want to perform Photospheric correction
    # --- check how to apply alphas
    if PhotosphericCorr:

        # --- if both alphas arent None then fall back to default (1,1)
        if alpha1 is not None and alpha2 is not None:
            alpha1 = 1
            alpha2 = 1

        # --- only alpha1 set -> get Pollux for alpha2
        elif alpha1 is not None and alpha2 is None:
            alpha2 = PolluxSpectra(
                        year=year,
                        filenum=filenum,
                        whichStar="Secondary",
                        whichSpectrum=alpha2_whichSpectrum,
                        fullResSpectra=False,
                        showHelp=False
                        )

        # --- only alpha2 set -> get Pollux for alpha1
        elif alpha2 is not None and alpha1 is None:
            alpha1 = PolluxSpectra(
                        year=year,
                        filenum=filenum,
                        whichStar="Primary",
                        whichSpectrum=alpha1_whichSpectrum,
                        fullResSpectra=False,
                        showHelp=False
                        )

        # --- both None -> Pollux for both
        else:
            alpha1 = PolluxSpectra(
                        year=year,
                        filenum=filenum,
                        whichStar="Primary",
                        whichSpectrum=alpha1_whichSpectrum,
                        fullResSpectra=False,
                        showHelp=False
                        )

            alpha2 = PolluxSpectra(
                        year=year,
                        filenum=filenum,
                        whichStar="Secondary",
                        whichSpectrum=alpha2_whichSpectrum,
                        fullResSpectra=False,
                        showHelp=False
                        )

    # --- if PhotosphericCorr is False then default alphas (1,1)
    else:
        alpha1 = 1
        alpha2 = 1
        beta1  = 0
        beta2  = 0

    # betas default
    if beta1 is None:
        beta1 = 0
    if beta2 is None:
        beta2 = 0

    # ---extract data
    LAM      = WAVE(year, filenum) # .................................................... lam
    F_LC     = NormFlux(year=year, filenum=filenum, returnFLC=True, showHelp=False) # ... F_LC_prime
    VTOT     = VISAMP(year, filenum, baseline) # ........................................ V_tot^Line
    VCONT    = VISAMPCONT(year, filenum, baseline) # .................................... V_tot^cont
    PHITOT   = np.deg2rad(VISPHITOT(year, filenum, baseline)) # ......................... PHI_TOT
    phitot   = np.deg2rad(VISPHI(year, filenum, baseline)) # ............................ phi_tot
    phicont  = np.deg2rad(DATA_CONT(year, filenum, baseline, dataType="VISPHI")) # ...... phi_cont
    U        = Extract(year, filenum, dataType="UCOORD")[baseline] # .................... U
    V        = Extract(year, filenum, dataType="VCOORD")[baseline] # .................... V
    DeltaX   = x_astrometrics[year_num[year]] * mas2rad # ............................... DeltaY
    DeltaY   = y_astrometrics[year_num[year]] * mas2rad # ............................... DeltaX
    phi2star = 2*np.pi / LAM * (U * DeltaX + V * DeltaY) # .............................. phi2star [rad]
    fc       = flux_ratios[year] # ...................................................... f_c

    # --- shortcut variables
    gamma     = (1 + beta1) * (1 + beta2) * (1 + fc)
    Gamma     = VTOT * F_LC
    D         = gamma * F_LC - (1 + beta2) * (alpha1 + beta1) - (1 + beta1) * (alpha2 + beta2) * fc
    S1        = (1 + beta2) * (1 - alpha1)
    S2        = (1 + beta1) * (1 - alpha2) * fc
    delta_phi = phicont - phi2star

    Term0 = gamma**2 * (Gamma**2 - 2 * Gamma * VCONT * np.cos(PHITOT) + VCONT**2)
    Term1 = 2 * gamma * S1 * (Gamma * np.cos(phitot) - VCONT * np.cos(phicont))
    Term2 = 2 * gamma * S2 * (Gamma * np.cos(phitot - phi2star) - VCONT * np.cos(phi2star - phicont))
    Term3 = S1**2 + S2**2 + 2 * S1 * S2 * np.cos(phi2star)

    # --- Pure line visibility under Photospheric absorption correction
    VIS_LINE_PA = 1/D * np.sqrt(Term0 + Term1 + Term2 + Term3)

    # --- compute error (old needs update)
    # --- define errors, we make some simplifications
    F_LC_ERR   = 0
    VCONT_ERR  = 0
    VTOT_ERR   = 0.01 * VTOT
    PHITOT_ERR = np.deg2rad(VISPHITOT(year, filenum, baseline)/VISPHITOT(year, filenum, baseline))

    term0 = 1/((F_LC - 1)**2 * VIS_LINE_PA)
    term1 = VCONT_ERR * (VCONT - F_LC * VTOT * np.cos(PHITOT))
    term2 = F_LC_ERR * (-F_LC * VTOT**2 + (F_LC +1) * VTOT * VCONT * np.cos(PHITOT) - VCONT**2)/(F_LC - 1)
    term3 = PHITOT_ERR * (F_LC * VTOT * VCONT * np.sin(PHITOT))
    term4 = VTOT_ERR * (F_LC**2 * VTOT - F_LC * VCONT * np.cos(PHITOT))

    VISLINE_ERR = term0 * np.sqrt(
        term1**2 + term2**2 + term3**2 + term4**2
    )

    # --- return error or Pure line visibility
    if returnError:
        return VISLINE_ERR
    else:
        return VIS_LINE_PA

# ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
# ┃                Mask for Brgamma interval           ┃
# ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

def mask_brg(year, filenum):

    # --- extract wavelength
    WL = WAVE(year, filenum)

    # --- define interval of bracket gamma emission
    left_interval  = np.array([2.1650, 2.1650, 2.1650, 2.1650, 2.1650, 2.1650, 2.1646])*1e-6
    right_interval = np.array([2.1679, 2.1696, 2.1679, 2.1682, 2.1682, 2.1682, 2.1680])*1e-6

    # --- define mask
    mask = (WL > left_interval[year_num[year]]) & (WL < right_interval[year_num[year]])

    return mask

# ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
# ┃              Pure Line Differential Phase          ┃
# ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
def PHILINE(
    year=None,
    filenum=None,
    baseline=None,
    wave=None,
    showHelp=False,
    plot=False,
    pltAllBaselines=False,
    PhotosphericCorr=False,
    alpha1_whichSpectrum=None,
    alpha2_whichSpectrum=None,
    alpha1=None,
    alpha2=None,
    beta1=None,
    beta2=None,
    returnError=False,
):

    # --- make sure that baseline is ignored if pltAllBaselines is True
    if pltAllBaselines:
        baseline = 0

    # --- Check if user inputed all necessary parameters
    if year is None or filenum is None or baseline is None:

        tbl = Table.grid(padding=(0, 2))
        tbl.add_column("Description", justify="right", style="bright_black")
        tbl.add_column("Usage", style=C_PARAM)

        tbl.add_row("choose Epoch......", f"[{C_PARAM}]year[/]=[{C_DEF}]epoch[/]")
        tbl.add_row("choose File number", f"[{C_PARAM}]filenum[/]=[{C_DEF}]value[/]")
        tbl.add_row("choose baseline...", f"[{C_PARAM}]baseline[/]=[{C_DEF}]value[/]")
        tbl.add_row("show help.........", f"[{C_PARAM}]showHelp[/]=[{C_DEF}]True[/]")

        panel = Panel(
            tbl,
            title="[bold #ff4747]ERROR --- Specify following parameters[/]",
            border_style=ERR_COLOR,
            box=box.ROUNDED,
            padding=(1, 2)
        )

        return console.print(panel)

    # --- Show help
    if showHelp:

        tbl = Table.grid(padding=(0, 2))
        tbl.add_column("Description", justify="right", style="bright_black")
        tbl.add_column("Usage", style=C_PARAM)

        tbl.add_row("apply Photospheric correction...........", f"[{C_PARAM}]PhotosphericCorr[/]=[{C_DEF}]True[/]")
        tbl.add_row("return Error............................", f"[{C_PARAM}]returnError[/]=[{C_DEF}]True[/]")
        tbl.add_row("choose which spectrum alpha1 uses.......", f"[{C_PARAM}]alpha1_whichSpectrum[/]=[{C_DEF}]value[/]")
        tbl.add_row("choose which spectrum alpha2 uses.......", f"[{C_PARAM}]alpha2_whichSpectrum[/]=[{C_DEF}]True[/]")
        tbl.add_row("choose alpha1 explicitely or set it to 1", f"[{C_PARAM}]alpha1[/]=[{C_DEF}]value[/]")
        tbl.add_row("choose alpha2 explicitely or set it to 1", f"[{C_PARAM}]alpha2[/]=[{C_DEF}]value[/]")
        tbl.add_row("choose beta1 explicitely or set it to 0.", f"[{C_PARAM}]beta1[/]=[{C_DEF}]value[/]")
        tbl.add_row("choose beta2 explicitely or set it to 0.", f"[{C_PARAM}]beta2[/]=[{C_DEF}]value[/]")
        tbl.add_row("plot data...............................", f"[{C_PARAM}]plot[/]=[{C_DEF}]True[/]")
        tbl.add_row("plot all baselines......................", f"[{C_PARAM}]pltAllBaselines[/]=[{C_DEF}]True[/]")
        tbl.add_row("change wavelength.......................", f"[{C_PARAM}]wave[/]=[{C_DEF}]value array[/]")

        panel = Panel(
            tbl,
            title="[bold #ffae42]Pure Line Differential Phase help[/]",
            border_style=C_FRAME,
            box=box.ROUNDED,
            padding=(1, 2)
        )

        console.print(panel)

    if pltAllBaselines:
        baseline2plot = [0, 1, 2, 3, 4, 5]
    else:
        baseline2plot = [baseline]

    for baseline in baseline2plot:

        if wave is not None:
            LAM = wave
        else:
            LAM = WAVE(year, filenum)

    # --- Check if user want to perform Photospheric correction
    # --- check how to apply alphas
    if PhotosphericCorr:

        # --- if both alphas arent None then fall back to default (1,1)
        if alpha1 is not None and alpha2 is not None:
            alpha1 = 1
            alpha2 = 1

        # --- only alpha1 set -> get Pollux for alpha2
        elif alpha1 is not None and alpha2 is None:
            alpha2 = PolluxSpectra(
                        year=year,
                        filenum=filenum,
                        whichStar="Secondary",
                        whichSpectrum=alpha2_whichSpectrum,
                        fullResSpectra=False,
                        returnSpectra=True,
                        showHelp=False
                        )

        # --- only alpha2 set -> get Pollux for alpha1
        elif alpha2 is not None and alpha1 is None:
            alpha1 = PolluxSpectra(
                        year=year,
                        filenum=filenum,
                        whichStar="Primary",
                        whichSpectrum=alpha1_whichSpectrum,
                        fullResSpectra=False,
                        returnSpectra=True,
                        showHelp=False
                        )

        # --- both None -> Pollux for both
        else:
            alpha1 = PolluxSpectra(
                        year=year,
                        filenum=filenum,
                        whichStar="Primary",
                        whichSpectrum=alpha1_whichSpectrum,
                        fullResSpectra=False,
                        returnSpectra=True,
                        showHelp=False
                        )

            alpha2 = PolluxSpectra(
                        year=year,
                        filenum=filenum,
                        whichStar="Secondary",
                        whichSpectrum=alpha2_whichSpectrum,
                        fullResSpectra=False,
                        returnSpectra=True,
                        showHelp=False
                        )

    # --- if PhotosphericCorr is False then default alphas (1,1)
    else:
        alpha1 = 1
        alpha2 = 1
        beta1  = 0
        beta2  = 0

    # betas default
    if beta1 is None:
        beta1 = 0
    if beta2 is None:
        beta2 = 0

    # --- extract data
    F_LC     = NormFlux(year=year, filenum=filenum, returnFLC=True) # ................ F_LC_prime
    VTOT     = VISAMP(year, filenum, baseline) # ..................................... V_tot^Line
    VCONT    = VISAMPCONT(year, filenum, baseline) # ................................. V_tot^cont
    PHITOT   = np.deg2rad(VISPHITOT(year, filenum, baseline)) # ...................... PHI_TOT
    phitot   = np.deg2rad(VISPHI(year, filenum, baseline)) # ......................... phi_tot
    phicont  = np.deg2rad(DATA_CONT(year, filenum, baseline, dataType="VISPHI")) # ... phi_cont
    U        = Extract(year, filenum, dataType="UCOORD")[baseline] # ................. U
    V        = Extract(year, filenum, dataType="VCOORD")[baseline] # ................. V
    DeltaX   = x_astrometrics[year_num[year]] * mas2rad # ............................ DeltaY
    DeltaY   = y_astrometrics[year_num[year]] * mas2rad # ............................ DeltaX
    phi2star = 2*np.pi / LAM * (U * DeltaX + V * DeltaY) # ........................... phi2star [rad]
    fc       = flux_ratios[year] # ................................................... f_c
    mask     = mask_brg(year, filenum) # ............................................. mask_brg

    # --- shortcut variables
    gamma     = (1 + beta1) * (1 + beta2) * (1 + fc)
    Gamma     = VTOT * F_LC
    D         = gamma * F_LC - (1 + beta2) * (alpha1 + beta1) - (1 + beta1) * (alpha2 + beta2) * fc
    S1        = (1 + beta2) * (1 - alpha1)
    S2        = (1 + beta1) * (1 - alpha2) * fc
    delta_phi = phicont - phi2star

    # ---- extract pure line visibility
    VLINE = VISLINE(
                year=year,
                filenum=filenum,
                baseline=baseline,
                PhotosphericCorr=False
            )

    # --- compute PHI_LINE
    term0 = 1/(D * VLINE)
    term1 = gamma * Gamma * np.sin(PHITOT)
    term2 = - S1 * np.sin(phicont)
    term3 = S2 * np.sin(phi2star - phicont)

    PHI_LINE = np.arcsin(
        term0 * (term1 + term2 + term3)
    )

    # --- compute PHI_LINE Error
    phi_tot_err = 0.0174533 # -------------------------------- 1 degree error in radian
    vis_tot_err = 0.01 * VISAMP(year, filenum, baseline) # --- 2% error on VisAmp
    vline_err   = VISLINE(
                    year=year,
                    filenum=filenum,
                    baseline=baseline,
                    PhotosphericCorr=False,
                    showHelp=False,
                    returnError=True
                    )
    f_lc_err     = 0

    err_term0 = np.tan(PHI_LINE)
    err_term1 = phi_tot_err / np.tan(PHITOT)
    err_term2 = vis_tot_err / VTOT
    err_term3 = vline_err / VLINE
    err_term4 = f_lc_err / (F_LC - 1)

    PHI_LINE_ERR = np.abs(err_term0 * np.sqrt(err_term1**2 + err_term2**2 + err_term3**2 + err_term4**2))

    # --- plot data
    if plot:
        fig, ax = plt.subplots(figsize=(6.3, 3.4))

        ax.plot(
            WL,
            VISPHITOT(year, filenum, baseline),
            c=baseline_colors(baseline)
        )

        ax.errorbar(
            WL[mask],
            np.rad2deg(PHI_LINE[mask]),
            yerr=np.rad2deg(PHI_LINE_ERR[mask]),
            fmt="s",
            c=baseline_colors(baseline),
            markeredgecolor="black",
            capsize=2,
            label=baseline_name(year, baseline)
        )

        ax.set_title(f"{baseline_name(year, baseline)}")
        ax.set_ylim(-70, 70)
        ax.set_xlim(2.16e-6, 2.174e-6)
        plt.legend()

    else:

        # --- return error or value
        if returnError:
            return PHI_LINE_ERR
        else:
            return PHI_LINE



# ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
# ┃                print radial velocity               ┃
# ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
def print_bary(year=None, filenum=None):

    if year is None or filenum is None:

        tbl = Table.grid(padding=(0, 2))
        tbl.add_column("Description", justify="right", style="bright_black")
        tbl.add_column("Usage", style=C_PARAM)

        tbl.add_row("choose Epoch......", f"[{C_PARAM}]year[/]=[{C_DEF}]epoch[/]")
        tbl.add_row("choose File number", f"[{C_PARAM}]filenum[/]=[{C_DEF}]value[/]")

        panel = Panel(
            tbl,
            title="[bold #ff4747]ERROR --- Specify following parameters[/]",
            border_style=ERR_COLOR,
            box=box.ROUNDED,
            padding=(1, 2)
        )

        return console.print(panel)

    # --- Barycentric correction
    rv_2018 = np.array([-20.887, -20.907, -20.917, -20.935, -20.948, -20.965, -20.977, -20.993])
    rv_2019 = np.array([-23.105, -23.127, -23.138, -23.158])
    rv_2020 = np.array([7.256, 7.232, 7.193, 7.181])
    rv_2021 = np.array([-13.099, -13.127, -13.141, -13.184, -13.211, -13.224, -13.252, -13.265, -13.291, -13.304, -13.328])
    rv_2022a = np.array([-12.350, -12.377, -12.391, -12.406])
    rv_2022b = np.array([-14.980, -15.006, -15.019, -15.032])
    rv_2023 = np.array([-8.434, -8.463, -8.492, -8.506])

    rv_epoch = {
        "2018": rv_2018,
        "2019": rv_2019,
        "2020": rv_2020,
        "2021": rv_2021,
        "2022a": rv_2022a,
        "2022b": rv_2022b,
        "2023": rv_2023
    }

    return rv_epoch[year[:-4]][filenum]



# ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
# ┃                scatter epoch point pos             ┃
# ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
def scatter_epoch_point(sim_base, date, ax=None, label=None, **scatter_kw):

    """Scatter the predicted on-sky position of the companion at a given calendar date."""
    if ax is None:
        ax = plt.gca()

    yr_to_s = 3.154e7
    pi_mas  = 1.866
    # convert calendar date → simulation time [yr] relative to ref_date
    ref_date = datetime(2018, 3, 6)
    t_yr = (date - ref_date).total_seconds() / yr_to_s

    # integrate a copy to avoid advancing your main sim
    sim1 = sim_base.copy()
    sim1.integrate(t_yr)
    p0, p1 = sim1.particles[0], sim1.particles[1]

    # relative position in AU
    dRA_AU  = p1.x - p0.x      # +East
    dDec_AU = p1.y - p0.y      # +North

    # convert to mas using your parallax factor (mas per AU)
    x_mas = dRA_AU  * pi_mas
    y_mas = dDec_AU * pi_mas

    # your orbit line uses plt.plot(y_orbit, x_orbit, ...), so match that order
    pt = ax.scatter(y_mas, x_mas,
                    facecolors='none', edgecolors='C2', linewidths=1.2, s=64,
                    zorder=12, **scatter_kw)
    if label:
        ax.annotate(label, xy=(y_mas, x_mas), xytext=(6, 6),
                    textcoords='offset points', fontsize=8)
    return pt


# ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
# ┃                   Gaussian Fit                     ┃
# ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
def FitGaussian(
    year,
    filenum,
    showHelp=False,
    printFitParameters=False,
    plotLine=None,
    do_plot=True,
    return_dict=False,
):
    """
    If return_dict=False (default):
        - behaves as before and returns only v_mid (km/s, rounded).
    If return_dict=True:
        - returns a dict with all fit and derived quantities.
    If do_plot=False:
        - skips the Matplotlib pop-up.
    """
    plt.rcParams.update({"xtick.top": False})

    WL       = WAVE(year, filenum)
    NORMFLUX = NormFlux(year, filenum, alpha1=1, alpha2=1,
                        returnFLCcorr=True, showHelp=False)  # raw FLCcorr
    lambda0  = 2.16612e-6  # m
    c        = 299792458   # m/s
    c_km_s   = 299792.458  # km/s
    Å_to_m   = 1e-10

    if showHelp:
        print("==========================================================================")
        print("Print fitted Parameters.\033[96mprintFitParameters=True\033[0m")
        print("plot arbitrary line.....\033[96mplotLine=Value\033[0m")
        print("==========================================================================")

    # --- functions for secondary x-axis ---
    def wavelength_to_velocity(w):
        return (w - lambda0) / lambda0 * c / 1000  # m → km/s

    def velocity_to_wavelength(v):
        return lambda0 * (1 + v * 1000 / c)

    def double_gaussian_with_offset(x, A1, mu1, s1, A2, mu2, s2, C):
        g1 = A1 * np.exp(-(x - mu1)**2 / (2 * s1**2))
        g2 = A2 * np.exp(-(x - mu2)**2 / (2 * s2**2))
        return g1 + g2 + C

    # --- window on the barycentric-corrected wavelength (meters) ---
    mask    = (WL >= 2.161e-6) & (WL <= 2.172e-6)
    xA_data = WL[mask] * 1e10      # convert m -> Å for the fit
    y_data  = NORMFLUX[mask]

    # --- good initials/bounds in Å ---
    initial_guess = [0.30, 21661.0, 2.5,  0.30, 21671.0, 2.5,  1.00]
    bounds = (
        [0.0, 21645.0, 0.5,  0.0, 21666.0, 0.5,  0.7],   # lower
        [2.0, 21668.0, 8.0,  2.0, 21685.0, 8.0,  1.3]    # upper
    )

    # --- fit in Å ---
    popt, pcov = curve_fit(
        double_gaussian_with_offset,
        xA_data, y_data,
        p0=initial_guess,
        bounds=bounds,
        maxfev=50000
    )

    A1, mu1, s1, A2, mu2, s2, C = popt

    # --- plot (Å) if requested --------------------------------------------
    if do_plot:
        xxA   = np.linspace(xA_data.min(), xA_data.max(), 1200)
        fit_y = double_gaussian_with_offset(xxA, *popt)

        fig, ax = plt.subplots(figsize=(6.3, 3.2))

        if plotLine is not None:
            if plotLine > 1e3:
                ax.axvline(x=plotLine * 1e-10, c="red", label="plotLine")
            else:
                ax.axvline(x=plotLine, c="red", label="plotLine")

        # --- add top velocity axis ---
        ax = plt.gca()
        secax = ax.secondary_xaxis(
            "top",
            functions=(wavelength_to_velocity, velocity_to_wavelength)
        )
        secax.set_xlabel("Velocity relative to Brγ [km/s]")

        ax.set_ylim(0.8, 1.7)
        ax.set_xlim(2.161e-6, 2.172e-6)
        ax.plot(WL, NORMFLUX, c="C1", lw=2, alpha=0.8, label="Data")
        ax.axvline(2.16612e-6, color="0.4", linestyle="--", linewidth=1, label=r"Br$\gamma$")
        ax.plot(xxA * 1e-10, fit_y, c="black", lw=1, label="double-G fit")
        ax.set_xlabel("Wavelength [m]")
        ax.set_ylabel("Normalized flux")
        plt.legend(frameon=False)
        plt.tight_layout()
        plt.show()

    plt.rcParams.update({"xtick.top": True})

    # ==================  UNCERTAINTIES / DERIVED QUANTITIES  ==================
    perr = np.sqrt(np.diag(pcov))  # 1σ formal errors (absolute_sigma=False)

    A1_err, mu1_err, s1_err, A2_err, mu2_err, s2_err, C_err = perr

    # --- derived in Å ---
    k_fwhm = 2.0 * np.sqrt(2.0 * np.log(2.0))   # FWHM = k * sigma

    FWHM1_A     = k_fwhm * s1
    FWHM2_A     = k_fwhm * s2
    FWHM1_A_err = k_fwhm * s1_err
    FWHM2_A_err = k_fwhm * s2_err

    # peak separation Δ = mu2 - mu1 (with covariance)
    cov     = pcov
    var_dmu = cov[4, 4] + cov[1, 1] - 2.0 * cov[4, 1]
    dmu_A   = mu2 - mu1
    dmu_A_err = np.sqrt(max(var_dmu, 0.0))

    # midpoint m = (mu1 + mu2)/2 (with covariance)
    var_mid   = 0.25 * (cov[1, 1] + cov[4, 4] + 2.0 * cov[1, 4])
    mid_A     = 0.5 * (mu1 + mu2)
    mid_A_err = np.sqrt(max(var_mid, 0.0))

    # --- convert to meters ---
    FWHM1_m, FWHM1_m_err = FWHM1_A * Å_to_m, FWHM1_A_err * Å_to_m
    FWHM2_m, FWHM2_m_err = FWHM2_A * Å_to_m, FWHM2_A_err * Å_to_m
    dmu_m,   dmu_m_err   = dmu_A * Å_to_m,   dmu_A_err * Å_to_m
    mid_m,   mid_m_err   = mid_A * Å_to_m,   mid_A_err * Å_to_m

    # ==================  RICH PRINTING (OPTIONAL)  ============================
    if printFitParameters:
        tbl = Table.grid(padding=(0, 2))
        tbl.add_column("Quantity", justify="right", style="bright_black")
        tbl.add_column("Usage", style=C_PARAM)

        # --- section: fit parameters (Å units for μ, s)
        tbl.add_row("[bold]Fit parameters[/] (Å units for centers/widths)", "")
        tbl.add_row("A1",  f"[{C_PARAM}]{A1:.4g}[/] ± [{C_DEF}]{A1_err:.2g}[/]")
        tbl.add_row("μ1",  f"[{C_PARAM}]{mu1:.4f} Å[/] ± [{C_DEF}]{mu1_err:.4f} Å[/]")
        tbl.add_row("s1",  f"[{C_PARAM}]{s1:.4f} Å[/] ± [{C_DEF}]{s1_err:.4f} Å[/]")
        tbl.add_row("A2",  f"[{C_PARAM}]{A2:.4g}[/] ± [{C_DEF}]{A2_err:.2g}[/]")
        tbl.add_row("μ2",  f"[{C_PARAM}]{mu2:.4f} Å[/] ± [{C_DEF}]{mu2_err:.4f} Å[/]")
        tbl.add_row("s2",  f"[{C_PARAM}]{s2:.4f} Å[/] ± [{C_DEF}]{s2_err:.4f} Å[/]")
        tbl.add_row("C",   f"[{C_PARAM}]{C:.5f}[/] ± [{C_DEF}]{C_err:.2g}[/]")

        # --- section: derived
        tbl.add_row("", "")
        tbl.add_row("[bold]Derived[/]", "")
        tbl.add_row("FWHM1",        f"[{C_PARAM}]{FWHM1_A:.4f} Å[/] ± [{C_DEF}]{FWHM1_A_err:.4f} Å[/]")
        tbl.add_row("FWHM2",        f"[{C_PARAM}]{FWHM2_A:.4f} Å[/] ± [{C_DEF}]{FWHM2_A_err:.4f} Å[/]")
        tbl.add_row("Peak sep Δ",   f"[{C_PARAM}]{dmu_A:.4f} Å[/] ± [{C_DEF}]{dmu_A_err:.4f} Å[/]")
        tbl.add_row("Midpoint (μ̄)", f"[{C_PARAM}]{mid_A:.4f} Å[/] ± [{C_DEF}]{mid_A_err:.4f} Å[/]")

        panel = Panel(
            tbl,
            title="[bold #ffae42]Double-Gaussian fit report[/]",
            border_style=C_FRAME,
            box=box.ROUNDED,
            padding=(1, 2),
        )
        console.print(panel)

    # ==================  CHI² + VELOCITIES + EW  =============================
    err_mid_f       = mid_A_err
    err_peak_sep_f  = dmu_A_err
    resid           = y_data - double_gaussian_with_offset(xA_data, *popt)
    dof             = max(1, len(y_data) - len(popt))
    chi2_red        = np.sum(resid**2) / dof
    cov_eff         = pcov * chi2_red

    # refined errors on Δμ, midpoint (from scaled cov)
    var_dmu_eff = cov_eff[4, 4] + cov_eff[1, 1] - 2.0 * cov_eff[4, 1]
    dmu_A_err   = np.sqrt(max(var_dmu_eff, 0.0))

    var_mid_eff = 0.25 * (cov_eff[1, 1] + cov_eff[4, 4] + 2.0 * cov_eff[1, 4])
    # mid_A_err = np.sqrt(max(var_mid_eff, 0.0))  # mid_A_err already used above

    dmu_m_err = dmu_A_err * Å_to_m
    mid_m_err = mid_A_err * Å_to_m

    delta_v_err = c_km_s * err_peak_sep_f / lambda0
    v_mid_err   = np.round(c_km_s * (err_mid_f * 1e-10) / lambda0, 2)

    delta_v = c_km_s * (mu2 * 1e-10 - mu1 * 1e-10) / lambda0
    v_mid   = c_km_s * (0.5 * (mu1 + mu2) * 1e-10 - lambda0) / lambda0

    # --- Equivalent widths ---------------------------------------------------
    k_gauss = np.sqrt(2 * np.pi)
    EW1_A   = -k_gauss * A1 * s1
    EW2_A   = -k_gauss * A2 * s2
    EW_total_A = EW1_A + EW2_A

    EW1_A_err = k_gauss * np.sqrt((s1 * A1_err)**2 + (A1 * s1_err)**2)
    EW2_A_err = k_gauss * np.sqrt((s2 * A2_err)**2 + (A2 * s2_err)**2)
    EW_total_A_err = np.sqrt(EW1_A_err**2 + EW2_A_err**2)

    EW1_m        = EW1_A * Å_to_m
    EW2_m        = EW2_A * Å_to_m
    EW_total_m   = EW_total_A * Å_to_m
    EW1_m_err    = EW1_A_err * Å_to_m
    EW2_m_err    = EW2_A_err * Å_to_m
    EW_total_m_err = EW_total_A_err * Å_to_m

    err_A1A2 = A1 / A2 * np.sqrt((A1_err / A1)**2 + (A2_err / A2)**2)

    if printFitParameters:
        # second rich panel: diagnostics
        tbl2 = Table.grid(padding=(0, 2))
        tbl2.add_column("Quantity", justify="right", style="bright_black")
        tbl2.add_column("Value", style=C_PARAM)

        tbl2.add_row("Δv...",       f"[{C_PARAM}]{delta_v:.3f}[/] ± [{C_DEF}]{delta_v_err:.2f}[/] km/s")
        tbl2.add_row("vmid...",     f"[{C_PARAM}]{v_mid:.3f}[/] ± [{C_DEF}]{v_mid_err:.2f}[/] km/s")
        tbl2.add_row("EW₁...",      f"[{C_PARAM}]{EW1_A:.3f}[/] ± [{C_DEF}]{EW1_A_err:.2f}[/] Å")
        tbl2.add_row("EW₂...",      f"[{C_PARAM}]{EW2_A:.3f}[/] ± [{C_DEF}]{EW2_A_err:.2f}[/] Å")
        tbl2.add_row("EW_total...", f"[{C_PARAM}]{EW_total_A:.3f}[/] ± [{C_DEF}]{EW_total_A_err:.2f}[/] Å")
        tbl2.add_row("V/R...",      f"[{C_PARAM}]{A1/A2:.3f}[/] ± [{C_DEF}]{err_A1A2:.2f}[/]")

        panel2 = Panel(
            tbl2,
            title="[bold #ffae42]Spectral line diagnostics[/]",
            border_style=C_FRAME,
            box=box.ROUNDED,
            padding=(1, 2),
        )
        console.print(panel2)

    # ==================  RETURN =============================================
    if return_dict:
        return {
            "A1": A1, "A1_err": A1_err,
            "A2": A2, "A2_err": A2_err,
            "mu1_A": mu1, "mu1_A_err": mu1_err,
            "mu2_A": mu2, "mu2_A_err": mu2_err,
            "s1_A": s1,  "s1_A_err": s1_err,
            "s2_A": s2,  "s2_A_err": s2_err,
            "C": C, "C_err": C_err,
            "FWHM1_A": FWHM1_A, "FWHM1_A_err": FWHM1_A_err,
            "FWHM2_A": FWHM2_A, "FWHM2_A_err": FWHM2_A_err,
            "FWHM1_m": FWHM1_m, "FWHM1_m_err": FWHM1_m_err,
            "FWHM2_m": FWHM2_m, "FWHM2_m_err": FWHM2_m_err,
            "dmu_A": dmu_A, "dmu_A_err": dmu_A_err,
            "mid_A": mid_A, "mid_A_err": mid_A_err,
            "dmu_m": dmu_m, "dmu_m_err": dmu_m_err,
            "mid_m": mid_m, "mid_m_err": mid_m_err,
            "delta_v": float(delta_v), "delta_v_err": float(delta_v_err*1e-10),
            "v_mid": float(v_mid), "v_mid_err": float(v_mid_err),
            "EW1_A": EW1_A, "EW1_A_err": EW1_A_err,
            "EW2_A": EW2_A, "EW2_A_err": EW2_A_err,
            "EW_total_A": EW_total_A, "EW_total_A_err": EW_total_A_err,
            "EW1_m": EW1_m, "EW1_m_err": EW1_m_err,
            "EW2_m": EW2_m, "EW2_m_err": EW2_m_err,
            "EW_total_m": EW_total_m, "EW_total_m_err": EW_total_m_err,
            "V_over_R": A1 / A2, "V_over_R_err": err_A1A2,
            "chi2_red": chi2_red,
        }

    # old behaviour: just give v_mid [km/s]
    return np.round(v_mid, 4)



# ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
# ┃             Annotate function for orbit            ┃
# ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
def annotate_func(text, x, y, dx, dy):
    return plt.annotate(
        text,
        xy=(x, y),
        xytext=(x - dx, y + dy),
        textcoords="data",
        fontsize=6,
        ha="center",
        va="center",
        arrowprops=dict(
            arrowstyle='-', lw=0.8, color='black'
        ),
        bbox=dict(
            boxstyle='round,pad=0.15',
            fc='white', alpha=0.6
        )
    )


# ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
# ┃                    Orbit plot                      ┃
# ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
def orbit(
    sma=None,
    ecc=None,
    inc=None,
    omega=None,
    Omega=None,
    tau=None,
    mtot=None,
    plx=None,
    labels=True,
    disableHelp=False,
    astrometrics=True,
    printParameters=False,
    plot=True,
    pltPredictedPoints=False,
    pltPrediction=None,
    pltPhotoCenter=None,
    returnPhotocenter=None,
    pltRV=False,
    pltRVxlim=None,
    returnRV=False,
    faceOn=False
):

    # --- Help
    if not disableHelp:

        tbl = Table.grid(padding=(0, 2))
        tbl.add_column("Description", justify="right", style="bright_black")
        tbl.add_column("Usage", style=C_PARAM)

        tbl.add_row("Disable this help...................", f"[{C_PARAM}]disableHelp[/]=[{C_DEF}]True[/]")
        tbl.add_row("Disable astrometric plotting........", f"[{C_PARAM}]astrometrics[/]=[{C_DEF}]False[/]")
        tbl.add_row("Disable labels......................", f"[{C_PARAM}]labels[/]=[{C_DEF}]False[/]")
        tbl.add_row("Print orbital parameters............", f"[{C_PARAM}]printParameters[/]=[{C_DEF}]True[/]")
        tbl.add_row("Disable Plot........................", f"[{C_PARAM}]plot[/]=[{C_DEF}]False[/]")
        tbl.add_row("Plot predicted positions............", f"[{C_PARAM}]pltPredictedPoints[/]=[{C_DEF}]True[/]")
        tbl.add_row("Plot specific epoch.................", f"[{C_PARAM}]pltPrediction[/]=[{C_DEF}]⟦year, month, day⟧[/]")
        tbl.add_row("Plot continuum photocenter..........", f"[{C_PARAM}]pltPhotoCenter[/]=[{C_DEF}]year[/]")
        tbl.add_row("Return photocenter coordinates......", f"[{C_PARAM}]returnPhotocenter[/]=[{C_DEF}]year[/]")
        tbl.add_row("Plot Radial Velocity................", f"[{C_PARAM}]pltRV[/]=[{C_DEF}]True[/]")
        tbl.add_row("Radial Velocity xlim................", f"[{C_PARAM}]pltRVxlim[/]=[{C_DEF}]⟦value1, value2⟧[/]")
        tbl.add_row("Return Radial Velocity all epochs...", f"[{C_PARAM}]returnRV[/]=[{C_DEF}]True[/]")
        tbl.add_row("Show orbit face on..................", f"[{C_PARAM}]faceOn[/]=[{C_DEF}]True[/]")

        panel = Panel(
            tbl,
            title="[bold #ffae42]Orbit Function Parameters[/]",
            border_style=C_FRAME,
            box=box.ROUNDED,
            padding=(1, 2)
        )

        console.print(panel)

    # --- Use best orbital parameters I know by default
    if sma is None:
        sma = 31.091
    if ecc is None:
        ecc = 0.990
    if inc is None:
        inc = 1.457
    if omega is None:
        omega = (1.910 + np.pi) % (2*np.pi)
    if Omega is None:
        Omega = (4.603 + np.pi) % (2*np.pi)
    if tau is None:
        tau=0.9894 # this one works best
        # tau = 0.971
    if mtot is None:
        mtot = 14.817
    if plx is None:
        plx = 1.866


    if faceOn:

        plt.figure(figsize=(6.4,3.2))
        x = np.linspace(-sma, sma, 1000)
        y1 = np.sqrt(1-ecc**2) * np.sqrt(-x**2 + sma**2)

        plt.title("Face on Orbit")
        plt.scatter([-sma*ecc],[0], marker="*", label="Primary", zorder=10, c="red")
        plt.plot(x,y1, c="b")
        plt.plot(x,-y1, label="Model Orbit", c="b")
        plt.legend(loc="upper right", frameon=True)
        return


    # --- print orbital parameters
    if printParameters:
        console.print(f"sma....{sma} [au]", markup=False)
        console.print(f"ecc......{ecc} [---]", markup=False)
        console.print(f"inc.....{np.rad2deg(inc):.4g} [deg]", markup=False)
        console.print(f"omega...{np.rad2deg(omega):.4g} [deg]", markup=False)
        console.print(f"Omega...{np.rad2deg(Omega):.4g} [deg]", markup=False)
        console.print(f"mtot...{mtot} [Msol]", markup=False)
        console.print(f"plx.....{plx} [mas]", markup=False)
        console.print(f"tau.....{tau} [---]", markup=False)
        console.print(f"P........{np.sqrt(sma**3/mtot):.4g} [yr]", markup=False)

    # --- return the coordinates of the photocenter
    if returnPhotocenter is not None:
        f_cont = flux_ratios[returnPhotocenter]
        PhotoCenter_x = f_cont/(1 + f_cont) * x_astrometrics[year_num[returnPhotocenter]]
        PhotoCenter_y = f_cont/(1 + f_cont) * y_astrometrics[year_num[returnPhotocenter]]

        return (PhotoCenter_x, PhotoCenter_y)


    # --- setup sim
    sim = rebound.Simulation()
    sim.units = ('yr', 'AU', 'Msun')

    m1 = 10.29
    sim.add(m=m1)
    sim.add(
        m=mtot-m1,
        a=sma,
        e=ecc,
        inc=inc,
        omega=omega,
        Omega=Omega,
        T= (1 - tau) * np.sqrt(sma**3/mtot)
    )

    sim.move_to_com()

    # --- Integrate one period
    t = np.linspace(-10, sim.particles[1].P, 60000)
    x_au         = np.empty_like(t)
    y_au         = np.empty_like(t)
    rv_primary   = np.empty_like(t)
    rv_secondary = np.empty_like(t)

    for i, ti in enumerate(t):
        sim.integrate(ti)
        p0, p1 = sim.particles[0], sim.particles[1]
        x_au[i] = p1.x - p0.x
        y_au[i] = p1.y - p0.y

        # --- in km/s thats why 4.74047
        rv_primary[i]   = -p0.vz * 4.74047
        rv_secondary[i] = -p1.vz * 4.74047

    # --- AU -> mas conversion
    x_mas = x_au * plx
    y_mas = y_au * plx

    # --- enable plot
    if plot:

        fig, ax = plt.subplots(figsize=(6.4,3.2))

        # --- labels of epochs
        if labels:
            annotate_func(text="06-03-2018", x=x_astrometrics[0], y=y_astrometrics[0], dx=1, dy=2)
            annotate_func(text="18-03-2019", x=x_astrometrics[1], y=y_astrometrics[1], dx=2.5, dy=1.5)
            annotate_func(text="23-12-2020", x=x_astrometrics[2], y=y_astrometrics[2], dx=1, dy=1.5)
            annotate_func(text="09-02-2021", x=x_astrometrics[3], y=y_astrometrics[3], dx=-1, dy=-1.5)
            annotate_func(text="07-02-2022", x=x_astrometrics[4], y=y_astrometrics[4], dx=1, dy=1.5)
            annotate_func(text="14-02-2022", x=x_astrometrics[5], y=y_astrometrics[5], dx=-1, dy=-1.5)
            annotate_func(text="29-01-2023", x=x_astrometrics[6], y=y_astrometrics[6], dx=0, dy=-2)


        ax.invert_xaxis()

        plt.axhline(y=0, c="black", ls=":", lw=0.8)
        plt.axvline(x=0, c="black", ls=":", lw=0.8)
        plt.scatter([0], [0], marker="+", c="black", zorder=10, label="Primary")

        plt.xlabel(r"$\delta x$ [mas] $\leftarrow$ East")
        plt.ylabel(r"$\delta y$ [mas] North $\rightarrow$")
        plt.title("Model Orbit with Astrometrics")
        plt.plot(
            y_mas,
            x_mas
        )

    # --- plot photocenter
    if pltPhotoCenter is not None:

        f_cont = flux_ratios[pltPhotoCenter]
        PhotoCenter_x = f_cont/(1 + f_cont) * x_astrometrics[year_num[pltPhotoCenter]]
        PhotoCenter_y = f_cont/(1 + f_cont) * y_astrometrics[year_num[pltPhotoCenter]]

        plt.scatter(
            PhotoCenter_x,
            PhotoCenter_y,
            marker="+",
            c="purple",
            label="Cont. Photocenter"
        )

    # --- plot the predicted points in the model
    if pltPredictedPoints:
        scatter_epoch_point(sim, datetime(2018, 3, 6), ax=plt.gca())
        scatter_epoch_point(sim, datetime(2019, 3, 18), ax=plt.gca())
        scatter_epoch_point(sim, datetime(2020, 12, 23), ax=plt.gca())
        scatter_epoch_point(sim, datetime(2021, 2, 9), ax=plt.gca())
        scatter_epoch_point(sim, datetime(2022, 2, 7), ax=plt.gca())
        scatter_epoch_point(sim, datetime(2022, 2, 14), ax=plt.gca())
        scatter_epoch_point(sim, datetime(2023, 1, 29), ax=plt.gca())

    if pltPrediction is None:
        pass
    else:
        scatter_epoch_point(
            sim,
            datetime(pltPrediction[0], pltPrediction[1], pltPrediction[2]),
            ax=plt.gca(),
            label=f"{pltPrediction[0]}-{pltPrediction[1]}-{pltPrediction[2]}"
        )

    if astrometrics:
        plt.errorbar(
            x_astrometrics,
            y_astrometrics,
            xerr = x_err_astrometrics,
            yerr = y_err_astrometrics,
            fmt="+",
            markeredgecolor="r",
            label="astrometrics",
            capsize=2
        )

    plt.legend(frameon=True, loc="upper right")


    # --- Returns the index of the time of the radial velocity
    def idx_rv(year, slectDate=None, changeRef=None):

        if changeRef is None:
            REF = datetime(2018, 3, 6)
        else:
            REF = changeRef

        epoch2date = {
            "2018_red":  datetime(2018, 3, 6),
            "2019_red":  datetime(2019, 3, 18),
            "2020_red":  datetime(2020, 12, 23),
            "2021_red":  datetime(2021, 2, 9),
            "2022a_red": datetime(2022, 2, 7),
            "2022b_red": datetime(2022, 2, 14),
            "2023_red":  datetime(2023, 1, 29)
        }

        if slectDate is None:

            # --- returns amount of years relative to REF date
            date2float = (epoch2date[year] - REF).total_seconds() / (365.25 * 86400)
            return np.argmin(np.abs(t - date2float))
        else:
            date2float = (slectDate - REF).total_seconds() / (365.25 * 86400)
            return np.argmin(np.abs(t - date2float))

    # --- print the radial velocity for these epochs
    if returnRV:

        for epoch in ["2018_red", "2019_red", "2020_red", "2021_red", "2022a_red", "2022b_red", "2023_red"]:

            idx = idx_rv(epoch)
            print(f"{epoch[:-4]} RV = {rv_primary[idx]:.8g} km/s")

    if pltRV:

        plt.close("all")

        fig, ax = plt.subplots(figsize=(6.4,3.2))

        if pltRVxlim is not None:
            ax.set_xlim(pltRVxlim[0], pltRVxlim[1])

        ax.plot(
            t,
            rv_primary,
            label="Primary"
        )

        ax.plot(
            t,
            rv_secondary,
            label="Secondary"
        )

        for epoch in ["2018_red", "2019_red", "2020_red", "2021_red", "2022a_red", "2022a_red", "2023_red"]:
            ax.scatter(
                t[idx_rv(year=epoch)],
                rv_primary[idx_rv(year=epoch)],
                marker="s",
                c="black",
                zorder=100
            )

        ax.set_xlabel("Time since reference epoch (2018-03-06) [yr]")
        ax.set_ylabel("Radial velocity [km/s]")
        ax.set_title("Radial Velocity for primary and secondary component")

        ax.legend()
        plt.show()

    if returnPhotocenter and pltPhotoCenter is not None:
        plt.close()
        return return_pc


# ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
# ┃                Theoretical PHI_LINE                ┃
# ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
def PHI_THEO(
    year=None,
    filenum=None,
    baseline=None,
    unit=None,
    lam=None
):

    # --- error log
    if year is None or filenum is None or baseline is None or unit is None:

        tbl = Table.grid(padding=(0, 2))
        tbl.add_column("Description", justify="right", style="bright_black")
        tbl.add_column("Usage", style=C_PARAM)

        tbl.add_row("choose Epoch..........", f"[{C_PARAM}]year[/]=[{C_DEF}]epoch[/]")
        tbl.add_row("choose File number....", f"[{C_PARAM}]filenum[/]=[{C_DEF}]value[/]")
        tbl.add_row("choose Baselin........", f"[{C_PARAM}]baseline[/]=[{C_DEF}]value[/]")
        tbl.add_row("choose unit rad or deg", f"[{C_PARAM}]unit[/]=[{C_DEF}]\"rad\" or \"deg\"[/]")

        panel = Panel(
            tbl,
            title="[bold #ff4747]ERROR --- Specify following parameters[/]",
            border_style=ERR_COLOR,
            box=box.ROUNDED,
            padding=(1, 2)
        )

        return console.print(panel)

    # --- extract U, V coordinates from baseline
    U = Extract(year, filenum, dataType="UCOORD")[baseline]
    V = Extract(year, filenum, dataType="VCOORD")[baseline]

    # obtain position of origin with respect to cont photo center
    xcoord, ycoord = orbit(returnPhotocenter=year, plot=False, disableHelp=True)
    P_x = -xcoord * mas2rad
    P_y = -ycoord * mas2rad

    if lam is None:
        mask = mask_brg(year, filenum)
        lam  = WAVE(year, filenum)[mask]

    PHI_THEO_NEW = - 2 * np.pi/lam * (U * P_x + V * P_y)

    if unit=="rad":
        return PHI_THEO_NEW
    elif unit=="deg":
        return np.rad2deg(PHI_THEO_NEW)

# ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
# ┃                     Apply k                        ┃
# ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
def apply_k(
    year=None,
    filenum=None,
    baseline=None,
    PhotosphericCorr=False,
    alpha1_whichSpectrum=None,
    alpha2_whichSpectrum=None,
    beta1=None,
    beta2=None,
    k_values=None,
    parity=None,
    showHelp=False
):



    # --- Show help
    if showHelp:

        tbl = Table.grid(padding=(0, 2))
        tbl.add_column("Description", justify="right", style="bright_black")
        tbl.add_column("Usage", style=C_PARAM)

        tbl.add_row("Apply Photospheric correction...........", f"[{C_PARAM}]PhotosphericCorr[/]=[{C_DEF}]True[/]")
        tbl.add_row("choose which spectrum alpha1 uses.......", f"[{C_PARAM}]alpha1_whichSpectrum[/]=[{C_DEF}]value[/]")
        tbl.add_row("choose which spectrum alpha2 uses.......", f"[{C_PARAM}]alpha2_whichSpectrum[/]=[{C_DEF}]True[/]")
        tbl.add_row("choose beta1 explicitely or set it to 0.", f"[{C_PARAM}]beta1[/]=[{C_DEF}]value[/]")
        tbl.add_row("choose beta2 explicitely or set it to 0.", f"[{C_PARAM}]beta2[/]=[{C_DEF}]value[/]")
        tbl.add_row("plot data...............................", f"[{C_PARAM}]plot[/]=[{C_DEF}]True[/]")
        tbl.add_row("plot all baselines......................", f"[{C_PARAM}]pltAllBaselines[/]=[{C_DEF}]True[/]")
        tbl.add_row("change wavelength.......................", f"[{C_PARAM}]wave[/]=[{C_DEF}]value array[/]")

        panel = Panel(
            tbl,
            title="[bold #ffae42]help[/]",
            border_style=C_FRAME,
            box=box.ROUNDED,
            padding=(1, 2)
        )

        console.print(panel)


    # --- error log
    if year is None or filenum is None or baseline is None or parity is None:

        tbl = Table.grid(padding=(0, 2))
        tbl.add_column("Description", justify="right", style="bright_black")
        tbl.add_column("Usage", style=C_PARAM)

        tbl.add_row("choose Epoch....................................", f"[{C_PARAM}]year[/]=[{C_DEF}]epoch[/]")
        tbl.add_row("choose File number..............................", f"[{C_PARAM}]filenum[/]=[{C_DEF}]value[/]")
        tbl.add_row("choose Baselin..................................", f"[{C_PARAM}]baseline[/]=[{C_DEF}]value[/]")
        tbl.add_row("choose parity = \"odd\", \"even\" or \"compare\"", f"[{C_PARAM}]unit[/]=[{C_DEF}]\"rad\" or \"deg\"[/]")

        panel = Panel(
            tbl,
            title="[bold #ff4747]ERROR --- Specify following parameters[/]",
            border_style=ERR_COLOR,
            box=box.ROUNDED,
            padding=(1, 2)
        )

        return console.print(panel)

    # --- make sure beta is defined
    if beta1 is None:
        beta1 = 0
    if beta2 is None:
        beta2 = 0

    # --- define the theoretical and the data corrected for photo absorption
    PHI_theo = PHI_THEO(year, filenum, baseline, unit="rad")
    PHI_data = PHILINE(
                year=year,
                filenum=filenum,
                baseline=baseline,
                PhotosphericCorr=PhotosphericCorr,
                alpha1_whichSpectrum=alpha1_whichSpectrum,
                alpha2_whichSpectrum=alpha2_whichSpectrum,
                beta1=beta1,
                beta2=beta2
                )[mask_brg(year, filenum)]

    # --- define both cases
    odd_case  = 1/np.pi * (PHI_theo + PHI_data)
    even_case = 1/np.pi * (PHI_theo - PHI_data)

    # --- check parity argument
    if parity == "odd":
        return odd_case
    elif parity == "even":
        return even_case
    elif parity == "compare":

        tbl = Table.grid(padding=(0, 2))
        tbl.add_column("Description", justify="right", style="bright_black")
        tbl.add_column("Usage", style=C_PARAM)

        tbl.add_row("odd case..........", f"{np.round(odd_case, 4)}")
        tbl.add_row("")
        tbl.add_row("even case.........", f"{np.round(even_case, 4)}")

        panel = Panel(
            tbl,
            title=f"[bold #89d470]k values | baseline {baseline} | year {year[:-4]} | file {filenum}[/]",
            border_style="#89d470",
            box=box.ROUNDED,
            padding=(1, 2)
        )

        return console.print(panel)


# ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
# ┃                  PHI_LINE_CORR                     ┃
# ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
def PHILINE_CORR(
    year=None,
    filenum=None,
    baseline=None,
    k_values=None,
    PhotosphericCorr=None,
    alpha1_whichSpectrum=None,
    alpha2_whichSpectrum=None,
    beta1=None,
    beta2=None,
    plot=None,
    showHelp=False,
    wave=None
):

    # --- Show help
    if showHelp:

        tbl = Table.grid(padding=(0, 2))
        tbl.add_column("Description", justify="right", style="bright_black")
        tbl.add_column("Usage", style=C_PARAM)

        tbl.add_row("Apply Photospheric correction...........", f"[{C_PARAM}]PhotosphericCorr[/]=[{C_DEF}]True[/]")
        tbl.add_row("choose which spectrum alpha1 uses.......", f"[{C_PARAM}]alpha1_whichSpectrum[/]=[{C_DEF}]value[/]")
        tbl.add_row("choose which spectrum alpha2 uses.......", f"[{C_PARAM}]alpha2_whichSpectrum[/]=[{C_DEF}]True[/]")
        tbl.add_row("choose beta1 explicitely or set it to 0.", f"[{C_PARAM}]beta1[/]=[{C_DEF}]value[/]")
        tbl.add_row("choose beta2 explicitely or set it to 0.", f"[{C_PARAM}]beta2[/]=[{C_DEF}]value[/]")
        tbl.add_row("select k_values manually................", f"[{C_PARAM}]k_values[/]=[{C_DEF}]np.array(⟦k1, ..., k6⟧)[/]")

        panel = Panel(
            tbl,
            title="[bold #ffae42]help[/]",
            border_style=C_FRAME,
            box=box.ROUNDED,
            padding=(1, 2)
        )

        console.print(panel)

    # --- error log
    if year is None or filenum is None or baseline is None:

        tbl = Table.grid(padding=(0, 2))
        tbl.add_column("Description", justify="right", style="bright_black")
        tbl.add_column("Usage", style=C_PARAM)

        tbl.add_row("choose Epoch....................................", f"[{C_PARAM}]year[/]=[{C_DEF}]epoch[/]")
        tbl.add_row("choose File number..............................", f"[{C_PARAM}]filenum[/]=[{C_DEF}]value[/]")
        tbl.add_row("choose Baseline.................................", f"[{C_PARAM}]baseline[/]=[{C_DEF}]value[/]")

        panel = Panel(
            tbl,
            title="[bold #ff4747]ERROR --- Specify following parameters[/]",
            border_style=ERR_COLOR,
            box=box.ROUNDED,
            padding=(1, 2)
        )

        return console.print(panel)

    Phi_line = PHILINE(
        year=year,
        filenum=filenum,
        baseline=baseline,
        PhotosphericCorr=PhotosphericCorr,
        alpha1_whichSpectrum=alpha1_whichSpectrum,
        alpha2_whichSpectrum=alpha2_whichSpectrum,
        showHelp=False,
        wave=wave
    )

    # --- default k value settings
    if k_values is None:

        if year == "2018_red":
            k_values = np.array([0, 1, 2, 1, 1, 1])

        elif year == "2019_red":
            k_values = np.array([0, 0, 0, 0, 0, 0])

        elif year == "2020_red":
            k_values = np.array([1, 1, 0, np.nan, -1, np.nan])

        elif year == "2021_red":
            k_values = np.array([1, 1, 1, 1, 1, 0])

        elif year == "2022a_red":
            k_values = np.array([2, 1, 0, -1, np.nan, -1])

        elif year == "2022b_red":
            k_values = np.array([1, 2, 2, 1, 1, 0])

        elif year == "2023_red":
            k_values = np.array([0, 1, 2, 1, 2, 1])


    # --- PhiLine corrected
    if np.isnan(k_values[baseline]):
        return np.full_like(Phi_line, np.nan)
    else:
        return (-1)**(int(k_values[baseline])) * Phi_line  + np.pi * k_values[baseline]
####

# ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
# ┃               pure line photo shift                ┃
# ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
def pureLinePhotoShift(
    year=None,
    filenum=None,
    baseline=None,
    k_values=None,
    wave_number=None,
    PhotosphericCorr=None,
    alpha1_whichSpectrum=None,
    alpha2_whichSpectrum=None,
    beta1=None,
    beta2=None,
    plot=None,
    plotPC=False,
    showHelp=False,
    returnProjectedP=False,
    doFit=False,
    returnPC=False,
    IncludePhi_theo=False,
    returnProjectedPtheo=False,
    return_wave=False,
    wave=None
):


    # --- Show help
    if showHelp:

        tbl = Table.grid(padding=(0, 2))
        tbl.add_column("Description", justify="right", style="bright_black")
        tbl.add_column("Usage", style=C_PARAM)

        tbl.add_row("Apply Photospheric correction...........", f"[{C_PARAM}]PhotosphericCorr[/]=[{C_DEF}]True[/]")
        tbl.add_row("choose which spectrum alpha1 uses.......", f"[{C_PARAM}]alpha1_whichSpectrum[/]=[{C_DEF}]value[/]")
        tbl.add_row("choose which spectrum alpha2 uses.......", f"[{C_PARAM}]alpha2_whichSpectrum[/]=[{C_DEF}]True[/]")

        tbl.add_row("choose beta1 explicitely or set it to 0.", f"[{C_PARAM}]beta1[/]=[{C_DEF}]value[/]")
        tbl.add_row("choose beta2 explicitely or set it to 0.", f"[{C_PARAM}]beta2[/]=[{C_DEF}]value[/]")

        tbl.add_row("choose wave number (wavelength).........", f"[{C_PARAM}]wave_number[/]=[{C_DEF}]value[/]")
        tbl.add_row("return wavelengths used.................", f"[{C_PARAM}]return_wave[/]=[{C_DEF}]True[/]")
        tbl.add_row("select k_values manually................", f"[{C_PARAM}]k_values[/]=[{C_DEF}]np.array(⟦k1, ..., k6⟧)[/]")
        tbl.add_row("plot....................................", f"[{C_PARAM}]plot[/]=[{C_DEF}]True[/]")

        panel = Panel(
            tbl,
            title="[bold #ffae42]help[/]",
            border_style=C_FRAME,
            box=box.ROUNDED,
            padding=(1, 2)
        )

        console.print(panel)

    if wave_number is None:
        return console.print("Select wave_number=value")


    # --- store for each baseline
    store_projectedP     = []
    store_projectedP_err = []
    store_projectedPtheo = []

    # --- some variables
    mask          = mask_brg(year, filenum)

    # --- In case I want to use a custom wavelength
    if wave is None:
        Lam       = WAVE(year, filenum)[mask][wave_number] # --------- masked wavelength
        lam2vel   = c_km_s * (WAVE(year, filenum)[mask] - bg)/bg
    else:
        Lam       = wave[mask][wave_number]
        lam2vel   = c_km_s * (wave[mask] - bg)/bg

    UCOORD        = Extract(year, filenum, dataType="UCOORD") # ------ U coordinate
    VCOORD        = Extract(year, filenum, dataType="VCOORD") # ------ V coordinate
    PA            = np.arctan2(UCOORD, VCOORD) # --------------------- Position angle in radian

    # --- I need this to know what color to give the colorbar
    # --- it returns it in km/s automatically
    if return_wave:
        return lam2vel

    for baseline in range(6):

        # --- declairing variables
        philine_corr     = PHILINE_CORR( # --------------------------- corrected pure line differential phase
                            year=year,
                            filenum=filenum,
                            baseline=baseline,
                            k_values=k_values,
                            PhotosphericCorr=PhotosphericCorr,
                            alpha1_whichSpectrum=alpha1_whichSpectrum,
                            alpha2_whichSpectrum=alpha2_whichSpectrum,
                            beta1=beta1,
                            beta2=beta2,
                            plot=None,
                            showHelp=False,
                            wave=wave
                        )

        philine_corr_err = PHILINE( # ---------------------------------- pure line differential phase error
                            year=year,
                            filenum=filenum,
                            baseline=baseline,
                            PhotosphericCorr=False,
                            returnError=True,
                            showHelp=False,
                            wave=wave
                        )

        B_i              = np.sqrt(UCOORD**2 + VCOORD**2)[baseline] # --- Baseline length

        # --- Projected B dot P = p_i * cos(PA - phi) in [mas]
        projectedP     = - philine_corr[mask][wave_number] / (2*np.pi) * Lam/B_i * rad2mas
        projectedP_err =   np.abs(philine_corr_err[mask][wave_number] / (2*np.pi) * Lam/B_i * rad2mas)
        store_projectedP.append(projectedP)
        store_projectedP_err.append(projectedP_err)

        # --- for the theoretical p_proj
        projectedPtheo = - PHI_THEO(
                            year,
                            filenum,
                            baseline,
                            unit="rad",
                        )[wave_number] / (2*np.pi) * Lam/B_i * rad2mas

        store_projectedPtheo.append(projectedPtheo)

    # --- transform into array since its a pain in the ass
    store_projectedP     = np.array(store_projectedP)
    store_projectedP_err = np.abs(np.array(store_projectedP_err))
    store_projectedPtheo = np.array(store_projectedPtheo)

    # return the projected photo shift
    if returnProjectedP:
        return store_projectedP

    if returnProjectedPtheo:
        return store_projectedPtheo

    # --- fit cosine angle and amplitude
    # --- define function that will be fitted
    def cosine(x, p, phi):
        return p * np.cos(x - phi)

    mask                = mask_brg(year, filenum) # --- mask interval
    x_data              = PA # ------------------------ radians
    y_data              = store_projectedP # ---------- mas
    y_data_err          = store_projectedP_err # ------ mas

    # --- baseline colors and names should match those of the remaining baselines
    set_baseline_labels = np.array(
        [baseline_name(year, baseline) for baseline in range(6)],
        dtype=object
    )
    set_baseline_colors = np.array(
        [baseline_colors(baseline) for baseline in range(6)],
        dtype=object
    )

    # --- fix issue with np.nan in data set
    nan_mask   = np.isfinite(x_data) & np.isfinite(y_data)
    x_data     = x_data[nan_mask]
    y_data     = y_data[nan_mask]
    y_data_err = y_data_err[nan_mask]
    set_baseline_labels = set_baseline_labels[nan_mask]
    set_baseline_colors = set_baseline_colors[nan_mask]

    # --- initial guess and bounds
    initial_guess = [3, np.deg2rad(50)]
    bounds = (
        [0.0, -2*np.pi], # -------- lower [p, phi]
        [8.0, 2*np.pi] # --- upper [p, phi]
    )

    # --- fit
    popt, pcov = curve_fit(
        cosine,
        x_data,
        y_data,
        sigma=y_data_err,
        absolute_sigma=False,
        p0=initial_guess,
        bounds=bounds,
        maxfev=50000
     )

    # --- fitted values
    p_fitted, phi_fitted = popt

    # --- error on fitted values
    perr = np.sqrt(np.diag(pcov))
    p_err, phi_err = perr

    # --- increase resolution on fitted cosine function
    x_highres = np.linspace(-2*np.pi, 2*np.pi, 100)
    y_highres = cosine(x_highres, p_fitted, phi_fitted)

    # --- Plot the fitted function and projected data
    if plot:

        #plt.figure(figsize=(6.3, 3.4))
        fig, ax = plt.subplots(figsize=(6.3,3.4))

        ax.set_xlabel("Position Angle [deg]")
        ax.set_ylabel("Proj. pure line pc. shift [mas]")
        ax.set_title(f"{year} File {filenum}")


        if IncludePhi_theo:

            # --- fit
            popt, pcov = curve_fit(
                cosine,
                PA,
                store_projectedPtheo,
                p0=initial_guess,
                bounds=bounds,
                maxfev=50000
             )

            # --- fitted values
            p_theo_fitted, phi_theo_fitted = popt

            # --- increase resolution on fitted cosine function
            x_highres      = np.linspace(-2*np.pi, 2*np.pi, 100)
            y_theo_highres = cosine(x_highres, p_theo_fitted, phi_theo_fitted)

            print(f'p_theo = {p_theo_fitted:.3f} ... phi_theo {np.rad2deg(phi_theo_fitted):.3f}')

            ax.plot(
                np.rad2deg(x_highres),
                y_theo_highres,
                c='red',
                ls="--",
                alpha=0.25,
                zorder=-1
            )

            ax.errorbar(
                np.rad2deg(PA),
                store_projectedPtheo,
                fmt="o",
                c="red",
                alpha=0.4
            )

            ax.errorbar(
                np.rad2deg(PA) + 180,
                -store_projectedPtheo,
                fmt="o",
                c="red",
                alpha=0.4
            )



        # --- stupid python loop to properly color the data
        for i in range(len(set_baseline_colors)):
            ax.errorbar(
                np.rad2deg(x_data[i]),#np.rad2deg(PA)[i],
                y_data[i],#store_projectedP[i],
                yerr=y_data_err[i],#store_projectedP_err[i],
                fmt="s",
                label=set_baseline_labels[i],
                markerfacecolor=set_baseline_colors[i],
                markeredgecolor="black",#set_baseline_colors[i],
                ecolor=set_baseline_colors[i],
                capsize=2
            )

            ax.errorbar(
                np.rad2deg(x_data[i]) + 180,#np.rad2deg(PA)[i] + 180,
                -y_data[i],#-store_projectedP[i],
                yerr=y_data_err[i],#store_projectedP_err[i],
                fmt="s",
                markerfacecolor=set_baseline_colors[i],
                markeredgecolor="black",#set_baseline_colors[i],
                ecolor=set_baseline_colors[i],
                capsize=2
            )

        plt.legend()

        if doFit:

            ax.set_ylim(-6.5, 5)
            param_text = (
                rf"$p = {p_fitted:.2f} \pm {p_err:.2f}\,\mathrm{{mas}}$" "\n"
                rf"$\varphi = {np.rad2deg(phi_fitted):.1f} \pm {np.rad2deg(phi_err):.1f}\,^\circ$"
            )

            ax.text(
                0.02, 0.98, param_text,
                transform=ax.transAxes,
                ha="left", va="top",
                bbox=dict(boxstyle="round", facecolor="white", alpha=0.8, edgecolor="black"),
            )

            ax.set_title(rf"Pure-line photocenter shift vs. position angle | ({year}, File {filenum})")

            ax.plot(
                np.rad2deg(x_highres),
                y_highres,
                c='blue',
                lw=1.2
            )

            ax.legend(ncol=2, loc="lower left")

        # --- actually plot the pure line
        # ---photocentershift with respect to the orbit


        px = p_fitted * np.sin(phi_fitted)
        py = p_fitted * np.cos(phi_fitted)

        px_err = np.sqrt(
            (np.sin(phi_fitted))**2 * p_err**2 +
            (p_fitted * np.cos(phi_fitted))**2 * phi_err**2
        )

        py_err = np.sqrt(
            (np.cos(phi_fitted))**2 * p_err**2 +
            (p_fitted*np.sin(phi_fitted))**2 * phi_err**2
        )

        # --- return PC values
        if returnPC:
            plt.close()
            return px, px_err, py, py_err

            ax.errorbar(
                px,
                py,
                xerr=px_err,
                yerr=py_err,
                capsize=2
            )

        # doesnt work
        if plotPC:
            return

            #fig, ax = plt.subplots(figsize=(6.3,3.4))

            #ax.set_xlabel("Position Angle [deg]")
            #ax.set_ylabel("Proj. pure line pc. shift [mas]")
            #ax.set_title(f"{year} File {filenum}")

            #ax.errorbar(
                #px,
                #py,
                #xerr=px_err,
                #yerr=py_err,
                #fmt="s",
                #capsize=2
            #)


# ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
# ┃                p_i = P*cos(delta_i)                ┃
# ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
def projected_p(year, filenum, baseline, lam=None, PHI=None):

    mask = mask_brg(year, filenum)

    if lam is None:
        lam  = WAVE(year, filenum)[mask]

    if PHI is None:
        PHI_chosen = PHI_THEO(year, filenum, baseline, lam=lam, unit="rad")
    else:
        PHI_chosen = PHI

    UCOORD = Extract(year, filenum, dataType="UCOORD")
    VCOORD = Extract(year, filenum, dataType="VCOORD")

    B_i = np.sqrt(UCOORD**2 + VCOORD**2)[baseline]

    return - PHI_chosen/B_i * lam/(2*np.pi) * rad2mas


# ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
# ┃                   Fit p and ohi                    ┃
# ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
def p_phi_fit(
    year,
    filenum,
    wave_num,
    k=None,
    keep=None,
    returnFit=False,
    ChoosePHIunwrap=None, # PHI_UNWARP in radians
    ChooseLam=None
):

    theo_p = {
        "2018_red": 3.43,
        "2019_red": 0.01,
        "2020_red": 2.39,
        "2021_red": 2.44,
        "2022a_red": 3.52,
        "2022b_red": 3.6,
        "2023_red": 3.68
    }

    theo_phi = {
        "2018_red": 1.32,
        "2019_red": 2.57,
        "2020_red": 0.82,
        "2021_red": 0.84,
        "2022a_red": 0.92,
        "2022b_red": 0.92,
        "2023_red": 0.96
    }

    if k is None:
        if year == "2018_red":
            k = [0., 1., 2., 1., 1., 1.]
            #[0., 1., 2., 1., 2., 1.] 2018 perfect match
        elif year == "2019_red":
            k = [0., 0., 0., 0., 0., 0.]
        elif year == "2020_red":
            k = [1., 1., 0., 0., -1., -1.]
            # k = [1, 1, 0, 0, -1, -1]# #best
        elif year == "2021_red":
            k = [ 0.,  1.,  2., 1., 2., 1.]
        elif year == "2022a_red":
            # k = [2., 1., 0., 0., 0., -1.]
            k = [2., 1., 0., -1., -1., -1.] # from my organized notebook
        elif year == "2022b_red":
            k = [1., 2., 2., 1., 1., 0.]
        elif year == "2023_red":
            k = [0., 1., 2., 1., 2., 1.]

    if keep is None:
        if year == "2018_red":
            keep = np.array([True, True, True, True, False, True])
        elif year == "2019_red":
            keep = np.array([True, True, True, True, True, True])
        elif year == "2020_red":
            keep = np.array([True, True, True, False, True, False])
            #keep = np.array([True, True, True, False, True, True]) #best
        elif year == "2021_red":
            keep = np.array([True, True, True, True, True, True])
        elif year == "2022a_red":
            # keep = np.array([True, True, True, False, False, True])
            keep = np.array([True, True, True, False, True, True])
        elif year == "2022b_red":
            keep = np.array([True, True, True, True, True, False])
        elif year == "2023_red":
            keep = np.array([True, True, True, False, False, True])

    mask = mask_brg(year, filenum)

    if ChoosePHIunwrap is not None:
        PHI_UNWARP = ChoosePHIunwrap
    else:

        PHI_UNWARP = np.array(
            [PHILINE(year, filenum, baseline, wave=ChooseLam)[0][mask][wave_num] * (-1)**(k[baseline]) + np.pi * k[baseline] for baseline in range(6)]
        )

    # PHI_UNWARP_ERR = np.array(
    #     [PHILINE(year, filenum, baseline, wave=ChooseLam)[1][mask][wave_num] * (-1)**(k[baseline]) + np.pi * k[baseline] for baseline in range(6)]
    # )
    PHI_UNWARP_ERR = np.array([PHILINE(year, filenum, baseline, wave=ChooseLam)[1][mask][wave_num] for baseline in range(6)])

    UCOORD = Extract(year, filenum, dataType="UCOORD")
    VCOORD = Extract(year, filenum, dataType="VCOORD")
    Bp     = np.sqrt(UCOORD**2 + VCOORD**2)

    # --- compute p_i = - Phi/2pi lam/B [mas]
    if ChooseLam is not None:
        lam = ChooseLam[mask][wave_num]
    else:
        lam = WAVE(year, filenum)[mask][wave_num]


    x = np.arctan2(UCOORD, VCOORD) # Position angle in radian

    if ChoosePHIunwrap is not None:
        y = - PHI_UNWARP[:, wave_num]/(2 * np.pi) * lam/Bp * rad2mas # p_i in mas
        # y = np.array([projected_p(year, filenum, baseline, lam=WAVE(year, filenum)[mask][wave_num], PHI=PHI_UNWARP[baseline]) for baseline in range(6)])
    else:
        y = - PHI_UNWARP/(2 * np.pi) * lam/Bp * rad2mas

    x = x[keep]
    y = y[keep]


    # --- error on p_i = P_i cos(phi_i - PA_i) | careful with NAN
    y_err = lam/(2 * np.pi * Bp) * PHI_UNWARP_ERR * rad2mas * 2 # added 2 to increase the error slightly
    y_err = y_err[keep]

    # --- fitting cosine function
    def cos_fit(x, p, phi):
        return p * np.cos(phi - x)

    xData = x
    yData = y
    initial_guess = [2, np.pi/2]
    bounds = (
            [0.0, -2*np.pi],  # lower
            [8.0, 2*np.pi]    # upper
        )

    popt, pcov = curve_fit(
        cos_fit,
        xData, yData,
        p0=initial_guess,
        bounds=bounds,
        sigma=y_err,
        absolute_sigma=False,
        maxfev=50000
    )

    sigma_p   = np.sqrt(pcov[0,0])
    sigma_phi = np.sqrt(pcov[1,1])
    p_fit, phi_fit =  popt

    # --- errors on ppc_shift vector
    sigma_x_err = np.sqrt((np.sin(phi_fit) * sigma_p)**2 + (p_fit * np.cos(phi_fit) * sigma_phi)**2)
    sigma_y_err = np.sqrt((np.cos(phi_fit) * sigma_p)**2 + (p_fit * np.sin(phi_fit) * sigma_phi)**2)

    # plot
    name_baseline = np.array([baseline_name(year, baseline) for baseline in range(6)])[keep]
    colors = np.array([baseline_colors(baseline) for baseline in range(6)])[keep]

    x_line = np.linspace(-2*np.pi, 2*np.pi, 1001)

    x_theo = np.arctan2(UCOORD, VCOORD)
    y_theo = np.array([projected_p(year, filenum, baseline)[wave_num] for baseline in range(6)])


    if returnFit:
        return p_fit, phi_fit, sigma_x_err, sigma_y_err

    else:

        plt.plot(
            x_line,
            cos_fit(x_line, p_fit, phi_fit)
        )

        # --- scatter plots
        for baseline in range(np.sum(keep)):

            plt.errorbar(
                x[baseline],
                y[baseline],
                yerr=y_err[baseline],
                c=colors[baseline],
                fmt="s"
                # label=baseline_name(year, baseline)
            )

            plt.errorbar(
                x[baseline]+np.pi,
                -y[baseline],
                yerr=y_err[baseline],
                c=colors[baseline],
                fmt="s"
                )

            plt.scatter(
                x_theo[baseline],
                y_theo[baseline],
                c=colors[baseline],
                label=baseline_name(year, baseline),
                marker="x",
                zorder=100
            )

            plt.scatter(
                x_theo[baseline] + np.pi,
                -y_theo[baseline],
                c=colors[baseline],
                label=baseline_name(year, baseline),
                marker="x",
                zorder=100
            )

        plt.legend(
            ncol=6,
            loc="upper center",
            bbox_to_anchor=(0.5, -0.15)   # center below, adjust -0.15 to move lower/higher
        )
        plt.xlabel("PA [rad]")
        plt.ylabel(r"$p \cdot \cos{(\varphi - PA)}$ [mas] ")
        plt.title(f"p_fit = {p_fit:.4g} pm {sigma_p:.3g} | phi_fit = {np.rad2deg(phi_fit):.4g} pm {np.rad2deg(sigma_phi):.3g} | wave_num {wave_num} | {year}")


# ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
# ┃                  brute k fitter                    ┃
# ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
def BruteKFitter(year, baseline, k_range=None):

    if k_range is None:
        k = [-3, -2, -1, 0, 1, 2, 3]

# ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
# ┃                     Plot vec(p)                    ┃
# ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
def plotPureLineShift(year, baseline, p, phi):
    return 0
