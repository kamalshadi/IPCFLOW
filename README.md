# IPC FLOW/PLOT

# Run the IPC flow GUI
### Step 1 - setup parameters
| Parameter | Description | Default |
| --- | --- | --- |
Segment | Segment length for IPC analysis | 768 |
Overlap | Overlap between data segments | 0.9 |
fcut | The largest frequency of interest | 50hz |

> IPC will be calculated with temporal resolution of $~\Delta t=\frac{(1-overlap)\times segment}{fs}$ and frequency resolution of $\Delta f =\frac{fs}{Segment}$. With the default setting we get $~\Delta t=15ms$$ and $~\Delta f=6Hz$

### Step 2 - setup postprocessing parameters  
__IPC Flow__ will collapse all the IPC value in the frequency band of interEst into a single time-series. This frequency band can be set in the **step 2** of the GUI. The time-series is further downsampled to the resoltion set by __timebin__ using __median__, __mean__, and __max__. The original time-series without downsampling is a dumped in __.mat__ format with suffix __-timebin__ while the downsampled data is dumped in __csv__ format for easy parsing with __-timebin__ suffix again.

The columns of the CSV file are:

1 . **Time bin index**  
2 . **Channel 1**  
3 . **Channel 2**  
4 . **Peak IPC**  
5 . **Mean IPC**  
6 . **Median IPC**  

> All these file are in the same study folder with the original filename and added suffix.

###### Special Note on Cyclic Rythm of IPC in The Band
---
IPC flow also calculates the rythmiticity of IPC in the frequency band of interest in three modes:
1. Pre-stimulus  
2. Post-stimulus  
3. Whole dataset  

The result of this analysis is logged as a __CSV__ file in the same study folder with suffix **-cyclic**

1. **mode**  
2. **Channel 1**  
3. **Channel 1**  
4. **Peak**  
5. **freq**  

### Step 3 - Input setting  

You can run IPC flow in two modes:  
1. Single file: Select the __vhdr__ file. Note if the **resting state** flag is on IPC flow expects to find 3 __vhdr__ files with **pre, mid, and post** suffix in the folder like:
```TE_C10_resting_pre.vhdr, TE_C10_resting_mid.vhdr, TE_C10_resting_post.vhdr```.

> Note that you can run Batch mode either with resting state data or TMS data depending on the rest flag status in the GUI.

2. Batch mode: Selct folders with all __vhdr__ files of interest. If resting flag is on, only files with **_resting** in their filename wil be parsed.
