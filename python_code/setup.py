import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
import math
import matplotlib.lines as mlines
from matplotlib.text import OffsetFrom
from tqdm import tqdm
import tevar_heatmap as th
from tqdm import tqdm_notebook as pbar
from scipy.stats import pearsonr
from scipy.stats import ttest_ind
from functools import reduce
import seaborn as sns