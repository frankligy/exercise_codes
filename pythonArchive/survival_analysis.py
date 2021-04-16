'''
Below is a tutorial for performming simple survival analysis in python,
using lifeline package
'''

# single category
from lifelines.datasets import load_waltons
df = load_waltons() # returns a Pandas DataFrame

T = df['T']
E = df['E']

from lifelines import KaplanMeierFitter
kmf = KaplanMeierFitter()
kmf.fit(T,E)
a = kmf.survival_function_
b = kmf.cumulative_density_
c = kmf.plot_survival_function(at_risk_counts=True)


# two categories
import matplotlib.pyplot as plt
fig,ax = plt.subplots()
kmf = KaplanMeierFitter()
T,E = [],[]
for name, grouped_df in df.groupby('group'):
    T.append(grouped_df['T'].values)
    E.append(grouped_df['E'].values)
    kmf.fit(grouped_df["T"], grouped_df["E"], label=name)
    kmf.plot_survival_function(ax=ax)
from lifelines.statistics import logrank_test
results = logrank_test(T[0], T[1], E[0], E[1])
ax.text(x=0.05,y=0.05,s='Log-rank test: p-value is {:.2e}'.format(results.p_value),weight='bold')


# cox regression analysis
from lifelines.datasets import load_regression_dataset
regression_dataset = load_regression_dataset() # a Pandas DataFrame
from lifelines import CoxPHFitter
# Using Cox Proportional Hazards model
cph = CoxPHFitter()
cph.fit(regression_dataset, 'T', event_col='E')
cph.print_summary()
cph.plot()





