#### UI imports
import streamlit as st
import plotly.graph_objs as go
import plotly.express as px

#### Dataset and data manipulation imports
import pandas as pd
import seaborn as sns
import numpy as np

#### Importing data
from get_values import *

st.title("Tune GUI")


st.markdown(
    """ 
    This is the **set** tune
    """
)
#  DataFrame
df = getValues()

set_df = df[df['type'] == 'set']

edited_df = st.data_editor(
    set_df,
    key="1",
    column_config={
    "type": None
    },
    use_container_width=True,
    num_rows="dynamic"
)

st.markdown(
    """ 
    Once we process the set tune, this is the **actual** tune
    """
)

actual_df = df[df['type'] == 'actual']

edited_df = st.data_editor(
    actual_df,
    key="2",
    column_config={
    "type": None
    },
    use_container_width=True,
    num_rows="dynamic"
)

################################### REPLACE CODE UNDER WITH TABLE

df = getValues()

# plotting set & actual tunes
fig = px.scatter(df,
    x="x", 
    y="y", 
    color="time", 
    symbol="type",
)
fig.update_layout(
    title="Set & Actual Tunes",
    xaxis_title='Qh',
    yaxis_title='Qv',
    legend=dict(x=0, y=1, traceorder='normal', orientation='h')
)

st.plotly_chart(fig)